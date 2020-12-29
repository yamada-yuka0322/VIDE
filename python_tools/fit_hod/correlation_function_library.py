#+
#   VIDE -- Void IDentification and Examination -- ./python_tools/fit_hod/correlation_function_library.py
#   Copyright (C) 2010-2014 Guilhem Lavaux
#   Copyright (C) 2011-2014 P. M. Sutter
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; version 2 of the License.
# 
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program; if not, write to the Free Software Foundation, Inc.,
#   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#+
#Version 1.1
#LATEST MODIFICATION: 15/05/2013
#This file contains the functions needed to compute:
#1)-the 2pt correlation function
#2)-the 2pt cross-correlation function
from mpi4py import MPI
import numpy as np
import scipy.weave as wv
import sys,os
import time

###### MPI DEFINITIONS ######
comm=MPI.COMM_WORLD
nprocs=comm.Get_size()
myrank=comm.Get_rank()





################################################################################
#This functions computes the TPCF (2pt correlation function) 
#from an N-body simulation. It takes into account boundary conditions
#VARIABLES:
#pos_g: array containing the positions of the galaxies
#pos_r: array containing the positions of the random particles catalogue
#BoxSize: Size of the Box. Units must be equal to those of pos_r/pos_g
#DD_action: compute number of galaxy pairs from data or read them---compute/read
#RR_action: compute number of random pairs from data or read them---compute/read
#DR_action: compute number of galaxy-random pairs or read them---compute/read
#DD_name: file name to write/read galaxy-galaxy pairs results
#RR_name: file name to write/read random-random pairs results
#DR_name: file name to write/read galaxy-random pairs results
#bins: number of bins to compute the 2pt correlation function
#Rmin: minimum radius to compute the 2pt correlation function
#Rmax: maximum radius to compute the 2pt correlation function
#USAGE: at the end of the file there is a example of how to use this function
def TPCF(pos_g,pos_r,BoxSize,DD_action,RR_action,DR_action,
         DD_name,RR_name,DR_name,bins,Rmin,Rmax,verbose=False):

    #dims determined requiring that no more 8 adyacent subboxes will be taken
    dims=int(BoxSize/Rmax)
    dims2=dims**2; dims3=dims**3

    ##### MASTER #####
    if myrank==0:

        #compute the indexes of the halo/subhalo/galaxy catalogue
        Ng=len(pos_g)*1.0; indexes_g=[]
        coord=np.floor(dims*pos_g/BoxSize).astype(np.int32)
        index=dims2*coord[:,0]+dims*coord[:,1]+coord[:,2]
        for i in range(dims3):
            ids=np.where(index==i)[0]
            indexes_g.append(ids)
        indexes_g=np.array(indexes_g)

        #compute the indexes of the random catalogue
        Nr=len(pos_r)*1.0; indexes_r=[]
        coord=np.floor(dims*pos_r/BoxSize).astype(np.int32)
        index=dims2*coord[:,0]+dims*coord[:,1]+coord[:,2]
        for i in range(dims3):
            ids=np.where(index==i)[0]
            indexes_r.append(ids)
        indexes_r=np.array(indexes_r)


        #compute galaxy-galaxy pairs: DD
        if DD_action=='compute':
            DD=DDR_pairs(bins,Rmin,Rmax,BoxSize,dims,indexes1=indexes_g,
                             indexes2=None,pos1=pos_g,pos2=None)
            if verbose:
                print(DD)
                print(np.sum(DD))
            #write results to a file
            write_results(DD_name,DD,bins,'radial')
        else:
            #read results from a file
            DD,bins_aux=read_results(DD_name,'radial')
            if bins_aux!=bins:
                print('Sizes are different!')
                sys.exit()

        #compute random-random pairs: RR
        if RR_action=='compute':
            RR=DDR_pairs(bins,Rmin,Rmax,BoxSize,dims,indexes1=indexes_r,
                             indexes2=None,pos1=pos_r,pos2=None)
            if verbose:
                print(RR)
                print(np.sum(RR))
            #write results to a file
            write_results(RR_name,RR,bins,'radial')
        else:
            #read results from a file
            RR,bins_aux=read_results(RR_name,'radial')
            if bins_aux!=bins:
                print('Sizes are different!')
                sys.exit()

        #compute galaxy-random pairs: DR
        if DR_action=='compute':
            DR=DDR_pairs(bins,Rmin,Rmax,BoxSize,dims,indexes1=indexes_g,
                         indexes2=indexes_r,pos1=pos_g,pos2=pos_r)
            if verbose:
                print(DR)
                print(np.sum(DR))
            #write results to a file
            write_results(DR_name,DR,bins,'radial')
        else:
            #read results from a file
            DR,bins_aux=read_results(DR_name,'radial')
            if bins_aux!=bins:
                print('Sizes are different!')
                sys.exit()


        #final procesing
        bins_histo=np.logspace(np.log10(Rmin),np.log10(Rmax),bins+1)
        middle=0.5*(bins_histo[:-1]+bins_histo[1:])
        DD*=1.0; RR*=1.0; DR*=1.0

        r,xi_r,error_xi_r=[],[],[]
        for i in range(bins):
            if (RR[i]>0.0): #avoid divisions by 0
                xi_aux,error_xi_aux=xi(DD[i],RR[i],DR[i],Ng,Nr)
                r.append(middle[i])
                xi_r.append(xi_aux)
                error_xi_r.append(error_xi_aux)

        r=np.array(r)
        xi_r=np.array(xi_r)
        error_xi_r=np.array(error_xi_r)
        
        return r,xi_r,error_xi_r



    ##### SLAVES #####
    else:
        if DD_action=='compute':
            DDR_pairs(bins,Rmin,Rmax,BoxSize,dims,
                      indexes1=None,indexes2=None,pos1=None,pos2=None)
        if RR_action=='compute':
            DDR_pairs(bins,Rmin,Rmax,BoxSize,dims,
                      indexes1=None,indexes2=None,pos1=None,pos2=None)           
        if DR_action=='compute':
            DDR_pairs(bins,Rmin,Rmax,BoxSize,dims,
                      indexes1=None,indexes2=None,pos1=None,pos2=None)
################################################################################

################################################################################
#This functions computes the TPCCF (2pt cross-correlation function) 
#from an N-body simulation. It takes into account boundary conditions
#VARIABLES:
#pos_g1: array containing the positions of the galaxies1
#pos_g2: array containing the positions of the galaxies2
#pos_r: array containing the positions of the random particles catalogue
#BoxSize: Size of the Box. Units must be equal to those of pos_r/pos_g1/pos_g2
#DD_action: compute number of galaxy pairs from data or read them---compute/read
#RR_action: compute number of random pairs from data or read them---compute/read
#DR_action: compute number of galaxy-random pairs or read them---compute/read
#DD_name: file name to write/read galaxy-galaxy pairs results
#RR_name: file name to write/read random-random pairs results
#DR_name: file name to write/read galaxy-random pairs results
#bins: number of bins to compute the 2pt correlation function
#Rmin: minimum radius to compute the 2pt correlation function
#Rmax: maximum radius to compute the 2pt correlation function
#USAGE: at the end of the file there is a example of how to use this function
def TPCCF(pos_g1,pos_g2,pos_r,BoxSize,
          D1D2_action,D1R_action,D2R_action,RR_action,
          D1D2_name,D1R_name,D2R_name,RR_name,
          bins,Rmin,Rmax,verbose=False):          
          

    #dims determined requiring that no more 8 adyacent subboxes will be taken
    dims=int(BoxSize/Rmax)
    dims2=dims**2; dims3=dims**3

    ##### MASTER #####
    if myrank==0:

        #compute the indexes of the halo1/subhalo1/galaxy1 catalogue
        Ng1=len(pos_g1)*1.0; indexes_g1=[]
        coord=np.floor(dims*pos_g1/BoxSize).astype(np.int32)
        index=dims2*coord[:,0]+dims*coord[:,1]+coord[:,2]
        for i in range(dims3):
            ids=np.where(index==i)[0]
            indexes_g1.append(ids)
        indexes_g1=np.array(indexes_g1)

        #compute the indexes of the halo2/subhalo2/galaxy2 catalogue
        Ng2=len(pos_g2)*1.0; indexes_g2=[]
        coord=np.floor(dims*pos_g2/BoxSize).astype(np.int32)
        index=dims2*coord[:,0]+dims*coord[:,1]+coord[:,2]
        for i in range(dims3):
            ids=np.where(index==i)[0]
            indexes_g2.append(ids)
        indexes_g2=np.array(indexes_g2)

        #compute the indexes of the random catalogue
        Nr=len(pos_r)*1.0; indexes_r=[]
        coord=np.floor(dims*pos_r/BoxSize).astype(np.int32)
        index=dims2*coord[:,0]+dims*coord[:,1]+coord[:,2]
        for i in range(dims3):
            ids=np.where(index==i)[0]
            indexes_r.append(ids)
        indexes_r=np.array(indexes_r)


        #compute galaxy1-galaxy2 pairs: D1D2
        if D1D2_action=='compute':
            D1D2=DDR_pairs(bins,Rmin,Rmax,BoxSize,dims,indexes1=indexes_g1,
                           indexes2=indexes_g2,pos1=pos_g1,pos2=pos_g2)
            if verbose:
                print(D1D2)
                print(np.sum(D1D2))
            #write results to a file
            write_results(D1D2_name,D1D2,bins,'radial')
        else:
            #read results from a file
            D1D2,bins_aux=read_results(D1D2_name,'radial')
            if bins_aux!=bins:
                print('Sizes are different!')
                sys.exit()

        #compute galaxy1-random pairs: D1R
        if D1R_action=='compute':
            D1R=DDR_pairs(bins,Rmin,Rmax,BoxSize,dims,indexes1=indexes_g1,
                          indexes2=indexes_r,pos1=pos_g1,pos2=pos_r)
            if verbose:
                print(D1R)
                print(np.sum(D1R))
            #write results to a file
            write_results(D1R_name,D1R,bins,'radial')
        else:
            #read results from a file
            D1R,bins_aux=read_results(D1R_name,'radial')
            if bins_aux!=bins:
                print('Sizes are different!')
                sys.exit()

        #compute galaxy2-random pairs: D2R
        if D2R_action=='compute':
            D2R=DDR_pairs(bins,Rmin,Rmax,BoxSize,dims,indexes1=indexes_g2,
                          indexes2=indexes_r,pos1=pos_g2,pos2=pos_r)
            if verbose:
                print(D2R)
                print(np.sum(D2R))
            #write results to a file
            write_results(D2R_name,D2R,bins,'radial')
        else:
            #read results from a file
            D2R,bins_aux=read_results(D2R_name,'radial')
            if bins_aux!=bins:
                print('Sizes are different!')
                sys.exit()

        #compute random-random pairs: RR
        if RR_action=='compute':
            RR=DDR_pairs(bins,Rmin,Rmax,BoxSize,dims,indexes1=indexes_r,
                         indexes2=None,pos1=pos_r,pos2=None)
            if verbose:
                print(RR)
                print(np.sum(RR))
            #write results to a file
            write_results(RR_name,RR,bins,'radial')
        else:
            #read results from a file
            RR,bins_aux=read_results(RR_name,'radial')
            if bins_aux!=bins:
                print('Sizes are different!')
                sys.exit()


        #final procesing
        bins_histo=np.logspace(np.log10(Rmin),np.log10(Rmax),bins+1)
        middle=0.5*(bins_histo[:-1]+bins_histo[1:])

        inside=np.where(RR>0)[0]
        D1D2=D1D2[inside]; D1R=D1R[inside]; D2R=D2R[inside]; RR=RR[inside]
        middle=middle[inside]

        D1D2n=D1D2*1.0/(Ng1*Ng2)
        D1Rn=D1R*1.0/(Ng1*Nr)
        D2Rn=D2R*1.0/(Ng2*Nr)
        RRn=RR*2.0/(Nr*(Nr-1.0))
        
        xi_r=D1D2n/RRn-D1Rn/RRn-D2Rn/RRn+1.0

        return middle,xi_r



    ##### SLAVES #####
    else:
        if D1D2_action=='compute':
            DDR_pairs(bins,Rmin,Rmax,BoxSize,dims,
                      indexes1=None,indexes2=None,pos1=None,pos2=None)
        if D1R_action=='compute':
            DDR_pairs(bins,Rmin,Rmax,BoxSize,dims,
                      indexes1=None,indexes2=None,pos1=None,pos2=None)
        if D2R_action=='compute':
            DDR_pairs(bins,Rmin,Rmax,BoxSize,dims,
                      indexes1=None,indexes2=None,pos1=None,pos2=None)
        if RR_action=='compute':
            DDR_pairs(bins,Rmin,Rmax,BoxSize,dims,
                      indexes1=None,indexes2=None,pos1=None,pos2=None)           
################################################################################


################################################################################
#This function is used to compute the DD file (the number of random pairs in a
#random catalogue) that it is need for massive computation of the 2pt 
#correlation function
#from an N-body simulation. It takes into account boundary conditions
#VARIABLES:
#pos_r: array containing the positions of the random particles catalogue
#BoxSize: Size of the Box. Units must be equal to those of pos_r/pos_g
#RR_name: file name to write/read random-random pairs results
#bins: number of bins to compute the 2pt correlation function
#Rmin: minimum radius to compute the 2pt correlation function
#Rmax: maximum radius to compute the 2pt correlation function
#USAGE: at the end of the file there is a example of how to use this function
def DD_file(pos_r,BoxSize,RR_name,bins,Rmin,Rmax):

    #dims determined requiring that no more 8 adyacent subboxes will be taken
    dims=int(BoxSize/Rmax)
    dims2=dims**2; dims3=dims**3

    ##### MASTER #####
    if myrank==0:

        #compute the indexes of the random catalogue
        Nr=len(pos_r)*1.0; indexes_r=[]
        coord=np.floor(dims*pos_r/BoxSize).astype(np.int32)
        index=dims2*coord[:,0]+dims*coord[:,1]+coord[:,2]
        for i in range(dims3):
            ids=np.where(index==i)[0]
            indexes_r.append(ids)
        indexes_r=np.array(indexes_r)

        #compute random-random pairs: RR
        RR=DDR_pairs(bins,Rmin,Rmax,BoxSize,dims,indexes1=indexes_r,
                     indexes2=None,pos1=pos_r,pos2=None)
        print(RR)
        print(np.sum(RR))
        #write results to a file
        write_results(RR_name,RR,bins,'radial')

    ##### SLAVES #####
    else:
        DDR_pairs(bins,Rmin,Rmax,BoxSize,dims,
                  indexes1=None,indexes2=None,pos1=None,pos2=None)           
################################################################################






################################################################################
####### COMPUTE THE NUMBER OF PAIRS IN A CATALOG ####### (x,y,z) very fast
################################################################################
def DDR_pairs(bins,Rmin,Rmax,BoxSize,dims,indexes1,indexes2,pos1,pos2):

    dims2=dims**2; dims3=dims**3

    #we put bins+1. The last bin is only for pairs separated by r=Rmax
    pairs=np.zeros(bins+1,dtype=np.int64) 

    ##### MASTER #####
    if myrank==0:
        #Master sends the indexes and particle positions to the slaves
        for i in range(1,nprocs):
            comm.send(pos1,dest=i,tag=6)
            comm.send(pos2,dest=i,tag=7)
            comm.send(indexes1,dest=i,tag=8)
            comm.send(indexes2,dest=i,tag=9)

        #Masters distributes the calculation among slaves
        for subbox in range(dims3):
            b=comm.recv(source=MPI.ANY_SOURCE,tag=1)
            comm.send(False,dest=b,tag=2)
            comm.send(subbox,dest=b,tag=3)

        #Master gathers partial results from slaves and returns the final result
        for j in range(1,nprocs):
            b=comm.recv(source=MPI.ANY_SOURCE,tag=1)
            comm.send(True,dest=b,tag=2)
            pairs_aux=comm.recv(source=b,tag=10)
            pairs+=pairs_aux

        #the last element is just for situations in which r=Rmax
        pairs[bins-1]+=pairs[bins]

        return pairs[:-1]


    ##### SLAVES #####
    else:
        #position of the center of each subbox
        sub_c=np.empty(3,dtype=np.float32)

        #slaves receive the positions and indexes
        pos1=comm.recv(source=0,tag=6)
        pos2=comm.recv(source=0,tag=7)
        indexes1=comm.recv(source=0,tag=8)
        indexes2=comm.recv(source=0,tag=9)

        comm.send(myrank,dest=0,tag=1)
        final=comm.recv(source=0,tag=2)
        while not(final):
            
            subbox=comm.recv(source=0,tag=3)
            core_ids=indexes1[subbox] #ids of the particles in the subbox
            pos0=pos1[core_ids]

            sub_c[0]=(subbox/dims2+0.5)*BoxSize/dims
            sub_c[1]=((subbox%dims2)/dims+0.5)*BoxSize/dims
            sub_c[2]=((subbox%dims2)%dims+0.5)*BoxSize/dims

            #galaxy-galaxy or random-random case
            if pos2==None: 
                #first: distances between particles in the same subbox
                distances_core(pos0,BoxSize,bins,Rmin,Rmax,pairs)

                #second: distances between particles in the subbox and particles around
                ids=indexes_subbox_neigh(sub_c,Rmax,dims,BoxSize,indexes1,subbox)
                if ids!=[]:
                    posN=pos1[ids]
                    DR_distances(pos0,posN,BoxSize,bins,Rmin,Rmax,pairs)

            #galaxy-random case
            else:          
                ids=indexes_subbox(sub_c,Rmax,dims,BoxSize,indexes2)
                posN=pos2[ids]
                DR_distances(pos0,posN,BoxSize,bins,Rmin,Rmax,pairs)

            comm.send(myrank,dest=0,tag=1)
            final=comm.recv(source=0,tag=2)

        print('cpu ',myrank,' finished: transfering data to master')
        comm.send(pairs,dest=0,tag=10)
################################################################################


################################################################################
#this function computes the distances between all the particles-pairs and
#return the number of pairs found in each distance bin
def distances_core(pos,BoxSize,bins,Rmin,Rmax,pairs):

    l=pos.shape[0]

    support = """
       #include <iostream>
       using namespace std;
    """
    code = """
       float middle=BoxSize/2.0;
       float dx,dy,dz,r;
       float x1,y1,z1,x2,y2,z2;
       float delta=log10(Rmax/Rmin)/bins;
       int bin,i,j;

       for (i=0;i<l;i++){
            x1=pos(i,0);
            y1=pos(i,1);
            z1=pos(i,2);
            for (j=i+1;j<l;j++){
                x2=pos(j,0);
                y2=pos(j,1);
                z2=pos(j,2);
                dx=(fabs(x1-x2)<middle) ? x1-x2 : BoxSize-fabs(x1-x2);
                dy=(fabs(y1-y2)<middle) ? y1-y2 : BoxSize-fabs(y1-y2);
                dz=(fabs(z1-z2)<middle) ? z1-z2 : BoxSize-fabs(z1-z2);
                r=sqrt(dx*dx+dy*dy+dz*dz);

               if (r>=Rmin && r<=Rmax){
                   bin=(int)(log10(r/Rmin)/delta);
                   pairs(bin)+=1; 
               }
            }   
       }
    """
    wv.inline(code,['pos','l','BoxSize','Rmin','Rmax','bins','pairs'],
              type_converters = wv.converters.blitz,
              support_code = support,libraries = ['m'])

    return pairs
################################################################################
#pos1---an array of positions
#pos2---an array of positions
#the function returns the number of pairs in distance bins between pos1 and pos2
def DR_distances(p1,p2,BoxSize,bins,Rmin,Rmax,pairs):

    l1=p1.shape[0]
    l2=p2.shape[0]

    support = """
       #include <iostream>
       using namespace std;
    """
    code = """
       float middle=BoxSize/2.0;
       float dx,dy,dz,r;
       float x1,y1,z1,x2,y2,z2;
       float delta=log10(Rmax/Rmin)/bins;
       int bin,i,j;

       for (i=0;i<l1;i++){
           x1=p1(i,0); 
           y1=p1(i,1);
           z1=p1(i,2);
           for (j=0;j<l2;j++){
               x2=p2(j,0); 
               y2=p2(j,1);
               z2=p2(j,2);
               dx=(fabs(x1-x2)<middle) ? x1-x2 : BoxSize-fabs(x1-x2);
               dy=(fabs(y1-y2)<middle) ? y1-y2 : BoxSize-fabs(y1-y2);
               dz=(fabs(z1-z2)<middle) ? z1-z2 : BoxSize-fabs(z1-z2);
               r=sqrt(dx*dx+dy*dy+dz*dz);

               if (r>=Rmin && r<=Rmax){
                   bin=(int)(log10(r/Rmin)/delta);
                   pairs(bin)+=1; 
               }
           }   
       }
    """
    wv.inline(code,['p1','p2','l1','l2','BoxSize','Rmin','Rmax','bins','pairs'],
              type_converters = wv.converters.blitz,
              support_code = support)

    return pairs
################################################################################


################################################################################
#this routine computes the IDs of all the particles within the neighboord cells
#that which can lie within the radius Rmax
def indexes_subbox(pos,Rmax,dims,BoxSize,indexes):

    #we add dims to avoid negative numbers. For example
    #if something hold between -1 and 5, the array to be
    #constructed should have indexes -1 0 1 2 3 4 5. 
    #To achieve this in a clever way we add dims
    i_min=int(np.floor((pos[0]-Rmax)*dims/BoxSize))+dims
    i_max=int(np.floor((pos[0]+Rmax)*dims/BoxSize))+dims
    j_min=int(np.floor((pos[1]-Rmax)*dims/BoxSize))+dims
    j_max=int(np.floor((pos[1]+Rmax)*dims/BoxSize))+dims
    k_min=int(np.floor((pos[2]-Rmax)*dims/BoxSize))+dims
    k_max=int(np.floor((pos[2]+Rmax)*dims/BoxSize))+dims
    
    i_array=np.arange(i_min,i_max+1)%dims
    j_array=np.arange(j_min,j_max+1)%dims
    k_array=np.arange(k_min,k_max+1)%dims

    PAR_indexes=np.array([])
    for i in i_array:
        for j in j_array:
            for k in k_array:
                num=dims**2*i+dims*j+k
                ids=indexes[num]
                PAR_indexes=np.concatenate((PAR_indexes,ids)).astype(np.int32)

    return PAR_indexes
################################################################################
#this routine returns the ids of the particles in the neighboord cells
#that havent been already selected
def indexes_subbox_neigh(pos,Rmax,dims,BoxSize,indexes,subbox):

    #we add dims to avoid negative numbers. For example
    #if something hold between -1 and 5, the array to be
    #constructed should have indexes -1 0 1 2 3 4 5. 
    #To achieve this in a clever way we add dims
    i_min=int(np.floor((pos[0]-Rmax)*dims/BoxSize))+dims
    i_max=int(np.floor((pos[0]+Rmax)*dims/BoxSize))+dims
    j_min=int(np.floor((pos[1]-Rmax)*dims/BoxSize))+dims
    j_max=int(np.floor((pos[1]+Rmax)*dims/BoxSize))+dims
    k_min=int(np.floor((pos[2]-Rmax)*dims/BoxSize))+dims
    k_max=int(np.floor((pos[2]+Rmax)*dims/BoxSize))+dims
    
    i_array=np.arange(i_min,i_max+1)%dims
    j_array=np.arange(j_min,j_max+1)%dims
    k_array=np.arange(k_min,k_max+1)%dims

    ids=np.array([])
    for i in i_array:
        for j in j_array:
            for k in k_array:
                num=dims**2*i+dims*j+k
                if num>subbox:
                    ids_subbox=indexes[num]
                    ids=np.concatenate((ids,ids_subbox)).astype(np.int32)
    return ids
################################################################################


################################################################################
#This function computes the correlation function and its error once the number
#of galaxy-galaxy, random-random & galaxy-random pairs are given together
#with the total number of galaxies and random points
def xi(GG,RR,GR,Ng,Nr):
    
    normGG=2.0/(Ng*(Ng-1.0))
    normRR=2.0/(Nr*(Nr-1.0))
    normGR=1.0/(Ng*Nr)

    GGn=GG*normGG
    RRn=RR*normRR
    GRn=GR*normGR
    
    xi=GGn/RRn-2.0*GRn/RRn+1.0

    fact=normRR/normGG*RR*(1.0+xi)+4.0/Ng*(normRR*RR/normGG*(1.0+xi))**2
    err=normGG/(normRR*RR)*np.sqrt(fact)
    err=err*np.sqrt(3.0)

    return xi,err
################################################################################









################################################################################
#This function writes partial results to a file
def write_results(fname,histogram,bins,case):
    f=open(fname,'w')
    if case=='par-perp':
        for i in range(len(histogram)):
            coord_perp=i/bins
            coord_par=i%bins
            f.write(str(coord_par)+' '+str(coord_perp)+' '+str(histogram[i])+'\n')
    elif case=='radial':
        for i in range(len(histogram)):
            f.write(str(i)+' '+str(histogram[i])+'\n')
    else:
        print('Error in the description of case:')
        print('Choose between: par-perp or radial')
    f.close()        
################################################################################
#This functions reads partial results of a file
def read_results(fname,case):

    histogram=[]

    if case=='par-perp':
        bins=np.around(np.sqrt(size)).astype(np.int64)

        if bins*bins!=size:
            print('Error finding the size of the matrix')
            sys.exit()

        f=open(fname,'r')
        for line in f.readlines():
            a=line.split()
            histogram.append(int(a[2]))
        f.close()
        histogram=np.array(histogram)
        return histogram,bins
    elif case=='radial':
        f=open(fname,'r')
        for line in f.readlines():
            a=line.split()
            histogram.append(int(a[1]))
        f.close()
        histogram=np.array(histogram)
        return histogram,histogram.shape[0]
    else:
        print('Error in the description of case:')
        print('Choose between: par-perp or radial')
################################################################################






############ EXAMPLE OF USAGE: TPCF ############
"""
points_g=150000
points_r=200000

BoxSize=500.0 #Mpc/h
Rmin=1.0      #Mpc/h
Rmax=50.0     #Mpc/h
bins=30

DD_action='compute'
RR_action='compute'
DR_action='compute'
DD_name='DD.dat'
RR_name='RR.dat'
DR_name='DR.dat'

if myrank==0:
    pos_g=np.random.random((points_g,3))*BoxSize
    pos_r=np.random.random((points_r,3))*BoxSize

    start=time.clock()
    r,xi_r,error_xi=TPCF(pos_g,pos_r,BoxSize,DD_action,RR_action,DR_action,
                         DD_name,RR_name,DR_name,bins,Rmin,Rmax,verbose=True)

    print r
    print xi_r
    print error_xi
    end=time.clock()
    print 'time:',end-start
else:
    pos_g=None; pos_r=None
    TPCF(pos_g,pos_r,BoxSize,DD_action,RR_action,DR_action,
         DD_name,RR_name,DR_name,bins,Rmin,Rmax,verbose=True)
"""


############ EXAMPLE OF USAGE: TPCCF ############
"""
points_g1=150000
points_g2=150000
points_r=200000

BoxSize=500.0 #Mpc/h
Rmin=1.0      #Mpc/h
Rmax=50.0     #Mpc/h
bins=30

D1D2_action='compute'; D1D2_name='D1D2.dat'
D1R_action='compute'; D1R_name='D1R.dat'
D2R_action='compute'; D2R_name='D2R.dat'
RR_action='compute'; RR_name='RR.dat'


if myrank==0:
    pos_g1=np.random.random((points_g1,3))*BoxSize
    pos_g2=np.random.random((points_g2,3))*BoxSize
    pos_r=np.random.random((points_r,3))*BoxSize

    r,xi_r=TPCCF(pos_g1,pos_g2,pos_r,BoxSize,
                 D1D2_action,D1R_action,D2R_action,RR_action,
                 D1D2_name,D1R_name,D2R_name,RR_name,
                 bins,Rmin,Rmax,verbose=True)

    print r
    print xi_r
else:
    pos_g1=None; pos_g2=None; pos_r=None
    TPCCF(pos_g1,pos_g2,pos_r,BoxSize,D1D2_action,D1R_action,D2R_action,
          RR_action,D1D2_name,D1R_name,D2R_name,RR_name,bins,Rmin,Rmax,
          verbose=True)     
"""

############ EXAMPLE OF USAGE: DD_file ############
"""
points_r=200000

BoxSize=500.0 #Mpc/h
Rmin=1.0      #Mpc/h
Rmax=50.0     #Mpc/h
bins=30

RR_name='RR.dat'

if myrank==0:
    pos_r=np.random.random((points_r,3))*BoxSize
    DD_file(pos_r,BoxSize,RR_name,bins,Rmin,Rmax)
else:
    pos_r=None
    DD_file(pos_r,BoxSize,RR_name,bins,Rmin,Rmax)
"""
