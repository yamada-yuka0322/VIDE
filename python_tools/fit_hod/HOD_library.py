#+
#   VIDE -- Void IDentification and Examination -- ./python_tools/fit_hod/HOD_library.py
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
import numpy as np
import readsnap
import readsubf
import sys
import time
import random

###############################################################################
#this function returns an array containing the positions of the galaxies (kpc/h)
#in the catalogue according to the fiducial density, M1 and alpha
#CDM halos with masses within [min_mass,max_mass], are populated
#with galaxies. The IDs and positions of the CDM particles belonging to the
#different groups are read from the snapshots
#If one needs to creates many catalogues, this function is not appropiate,
#since it wastes a lot of time reading the snapshots and sorting the IDs
#min_mass and max_mass are in units of Msun/h, not 1e10 Msun/h
#mass_criteria: definition of the halo virial radius -- 't200' 'm200' 'c200'
#fiducial_density: galaxy number density to be reproduced, in (h/Mpc)^3
def hod(snapshot_fname,groups_fname,groups_number,min_mass,max_mass,
        fiducial_density,M1,alpha,mass_criteria,verbose=False):

    thres=1e-3 #controls the max relative error to accept a galaxy density
    
    #read the header and obtain the boxsize
    head=readsnap.snapshot_header(snapshot_fname)
    BoxSize=head.boxsize    #BoxSize in kpc/h

    #read positions and IDs of DM particles: sort the IDs array
    DM_pos=readsnap.read_block(snapshot_fname,"POS ",parttype=-1) #kpc/h
    DM_ids=readsnap.read_block(snapshot_fname,"ID  ",parttype=-1)-1 
    sorted_ids=DM_ids.argsort(axis=0)
    #the particle whose ID is N is located in the position sorted_ids[N]
    #i.e. DM_ids[sorted_ids[N]]=N
    #the position of the particle whose ID is N would be:
    #DM_pos[sorted_ids[N]]

    #read the IDs of the particles belonging to the CDM halos
    halos_ID=readsubf.subf_ids(groups_fname,groups_number,0,0,
                               long_ids=True,read_all=True)
    IDs=halos_ID.SubIDs-1
    del halos_ID

    #read CDM halos information
    halos=readsubf.subfind_catalog(groups_fname,groups_number,
                                   group_veldisp=True,masstab=True,
                                   long_ids=True,swap=False)
    if mass_criteria=='t200':
        halos_mass=halos.group_m_tophat200*1e10   #masses in Msun/h
        halos_radius=halos.group_r_tophat200      #radius in kpc/h
    elif mass_criteria=='m200':
        halos_mass=halos.group_m_mean200*1e10     #masses in Msun/h
        halos_radius=halos.group_r_mean200        #radius in kpc/h
    elif mass_criteria=='c200':    
        halos_mass=halos.group_m_crit200*1e10     #masses in Msun/h
        halos_radius=halos.group_r_crit200        #radius in kpc/h
    else:
        print('bad mass_criteria')
        sys.exit()
    halos_pos=halos.group_pos   #positions in kpc/h
    halos_len=halos.group_len
    halos_offset=halos.group_offset
    halos_indexes=np.where((halos_mass>min_mass) & (halos_mass<max_mass))[0]
    del halos
    
    if verbose:
        print(' ')
        print('total halos found=',halos_pos.shape[0])
        print('halos number density=',len(halos_pos)/(BoxSize*1e-3)**3)

    #keep only the halos in the given mass range 
    halo_mass=halos_mass[halos_indexes]
    halo_pos=halos_pos[halos_indexes]
    halo_radius=halos_radius[halos_indexes]
    halo_len=halos_len[halos_indexes]
    halo_offset=halos_offset[halos_indexes]
    del halos_indexes

    ##### COMPUTE Mmin GIVEN M1 & alpha #####
    i=0; max_iterations=20 #maximum number of iterations
    Mmin1=min_mass; Mmin2=max_mass
    while (i<max_iterations):
        Mmin=0.5*(Mmin1+Mmin2) #estimation of the HOD parameter Mmin

        total_galaxies=0
        inside=np.where(halo_mass>Mmin)[0] #take all galaxies with M>Mmin
        mass=halo_mass[inside] #only halos with M>Mmin have central/satellites

        total_galaxies=mass.shape[0]+np.sum((mass/M1)**alpha)
        mean_density=total_galaxies*1.0/(BoxSize*1e-3)**3 #galaxies/(Mpc/h)^3

        if (np.absolute((mean_density-fiducial_density)/fiducial_density)<thres):
            i=max_iterations
        elif (mean_density>fiducial_density):
            Mmin1=Mmin
        else:
            Mmin2=Mmin
        i+=1

    if verbose:
        print(' ')
        print('Mmin=',Mmin)
        print('average number of galaxies=',total_galaxies)
        print('average galaxy density=',mean_density)
    #########################################

    #just halos with M>Mmin; the rest do not host central/satellite galaxies
    inside=np.where(halo_mass>Mmin)[0]
    halo_mass=halo_mass[inside]
    halo_pos=halo_pos[inside]
    halo_radius=halo_radius[inside]
    halo_len=halo_len[inside]
    halo_offset=halo_offset[inside]
    del inside

    #compute number of satellites in each halo using the Poisson distribution 
    N_mean_sat=(halo_mass/M1)**alpha #mean number of satellites
    N_sat=np.empty(len(N_mean_sat),dtype=np.int32)
    for i in range(len(N_sat)):
        N_sat[i]=np.random.poisson(N_mean_sat[i])
    N_tot=np.sum(N_sat)+len(halo_mass) #total number of galaxies in the catalogue

    if verbose:
        print(' ')
        print(np.min(halo_mass),'< M_halo <',np.max(halo_mass))
        print('total number of galaxies=',N_tot)
        print('galaxy number density=',N_tot/(BoxSize*1e-3)**3)

    #put satellites following the distribution of dark matter in groups
    if verbose:
        print(' ')
        print('Creating mock catalogue ...', end=' ')

    pos_galaxies=np.empty((N_tot,3),dtype=np.float32)
    #index: variable that go through halos (may be several galaxies in a halo)
    #i: variable that go through all (central/satellites) galaxies
    #count: find number of galaxies that lie beyond its host halo virial radius
    index=0; count=0; i=0 
    while (index<halo_mass.shape[0]):

        position=halo_pos[index]  #position of the DM halo
        radius=halo_radius[index] #radius of the DM halo

        #save the position of the central galaxy
        pos_galaxies[i]=position; i+=1

        #if halo contains satellites, save their positions
        Nsat=N_sat[index]
        if Nsat>0:
            offset=halo_offset[index] 
            length=halo_len[index]
            idss=sorted_ids[IDs[offset:offset+length]]

            #compute the distances to the halo center keeping those with R<Rvir
            pos=DM_pos[idss] #positions of the particles belonging to the halo
            posc=pos-position

            #this is to populate correctly halos closer to box boundaries
            if np.any((position+radius>BoxSize) + (position-radius<0.0)):
                
                inside=np.where(posc[:,0]>BoxSize/2.0)[0]
                posc[inside,0]-=BoxSize
                inside=np.where(posc[:,0]<-BoxSize/2.0)[0]
                posc[inside,0]+=BoxSize

                inside=np.where(posc[:,1]>BoxSize/2.0)[0]
                posc[inside,1]-=BoxSize
                inside=np.where(posc[:,1]<-BoxSize/2.0)[0]
                posc[inside,1]+=BoxSize

                inside=np.where(posc[:,2]>BoxSize/2.0)[0]
                posc[inside,2]-=BoxSize
                inside=np.where(posc[:,2]<-BoxSize/2.0)[0]
                posc[inside,2]+=BoxSize
                
            radii=np.sqrt(posc[:,0]**2+posc[:,1]**2+posc[:,2]**2)
            inside=np.where(radii<radius)[0]
            selected=random.sample(inside,Nsat)
            pos=pos[selected]

            #aditional, not esential check. Can be comment out
            posc=pos-position
            if np.any((posc>BoxSize/2.0) + (posc<-BoxSize/2.0)):
                inside=np.where(posc[:,0]>BoxSize/2.0)[0]
                posc[inside,0]-=BoxSize
                inside=np.where(posc[:,0]<-BoxSize/2.0)[0]
                posc[inside,0]+=BoxSize

                inside=np.where(posc[:,1]>BoxSize/2.0)[0]
                posc[inside,1]-=BoxSize
                inside=np.where(posc[:,1]<-BoxSize/2.0)[0]
                posc[inside,1]+=BoxSize

                inside=np.where(posc[:,2]>BoxSize/2.0)[0]
                posc[inside,2]-=BoxSize
                inside=np.where(posc[:,2]<-BoxSize/2.0)[0]
                posc[inside,2]+=BoxSize
            r_max=np.max(np.sqrt(posc[:,0]**2+posc[:,1]**2+posc[:,2]**2))
            if r_max>radius: #check no particles beyond Rv selected 
                print(position)
                print(radius)
                print(pos)
                count+=1

            for j in range(Nsat):
                pos_galaxies[i]=pos[j]; i+=1
        index+=1

    if verbose:
        print('done')
    #some final checks
    if i!=N_tot:
        print('some galaxies missing:')
        print('register',i,'galaxies out of',N_tot)
    if count>0:
        print('error:',count,'particles beyond the virial radius selected')

    return pos_galaxies
###############################################################################





#This function is equal to the above one, except that the snapshot read, halos
#read and ID sorting it is not performing here. It is best suited when many 
#galaxy catalogues need to be created: for example, when iterating among M1 and 
#alpha trying to find the best combination that reproduces the measured wp(r)
#VARIABLES:
#DM_pos: array containing the positions of the CDM particles
#sorted_ids: array containing the positions of the IDs in the snapshots. 
#sorted_ids[N] gives the position where the particle whose ID is N is located
#IDs:IDs array as read from the subfind ID file
#halo_mass: array containing the masses of the CDM halos in the mass interval
#halo_pos:  array containing the positions of the CDM halos in the mass interval
#halo_radius: array containing the radii of the CDM halos in the mass interval
#halo_len: array containing the len of the CDM halos in the mass interval
#halo_offset: array containing the offset of the CDM halos in the mass interval
#BoxSize: Size of the simulation Box. In Mpc/h
#fiducial_density: galaxy number density to be reproduced, in (h/Mpc)^3
def hod_fast(DM_pos,sorted_ids,IDs,halo_mass,halo_pos,halo_radius,halo_len,
             halo_offset,BoxSize,min_mass,max_mass,fiducial_density,
             M1,alpha,seed,verbose=False):

    problematic_cases=0 #number of problematic cases (e.g. halos with Rvir=0.0)
    thres=1e-3 #controls the max relative error to accept a galaxy density

    ##### COMPUTE Mmin GIVEN M1 & alpha #####
    i=0; max_iterations=20 #maximum number of iterations
    Mmin1=min_mass; Mmin2=max_mass
    while (i<max_iterations):
        Mmin=0.5*(Mmin1+Mmin2) #estimation of the HOD parameter Mmin

        total_galaxies=0
        inside=np.where(halo_mass>Mmin)[0]
        mass=halo_mass[inside] #only halos with M>Mmin have central/satellites

        total_galaxies=mass.shape[0]+np.sum((mass/M1)**alpha)
        mean_density=total_galaxies*1.0/BoxSize**3

        if (np.absolute((mean_density-fiducial_density)/fiducial_density)<thres):
            i=max_iterations
        elif (mean_density>fiducial_density):
            Mmin1=Mmin
        else:
            Mmin2=Mmin
        i+=1

    if verbose:
        print(' ')
        print('Mmin=',Mmin)
        print('average number of galaxies=',total_galaxies)
        print('average galaxy density=',mean_density)
    #########################################

    #just halos with M>Mmin; the rest do not host central/satellite galaxies
    inside=np.where(halo_mass>Mmin)[0]
    halo_mass=halo_mass[inside]
    halo_pos=halo_pos[inside]
    halo_radius=halo_radius[inside]
    halo_len=halo_len[inside]
    halo_offset=halo_offset[inside]
    del inside

    #compute number of satellites in each halo using the Poisson distribution 
    np.random.seed(seed) #this is just to check convergence on w_p(r_p)
    N_mean_sat=(halo_mass/M1)**alpha #mean number of satellites
    N_sat=np.empty(len(N_mean_sat),dtype=np.int32)
    for i in range(len(N_sat)):
        N_sat[i]=np.random.poisson(N_mean_sat[i])
    N_tot=np.sum(N_sat)+len(halo_mass) #total number of galaxies in the catalogue

    if verbose:
        print(' ')
        print(np.min(halo_mass),'< M_halo <',np.max(halo_mass))
        print('total number of galaxies=',N_tot)
        print('galaxy number density=',N_tot/BoxSize**3)

    #put satellites following the distribution of dark matter in groups
    if verbose:
        print(' ')
        print('Creating mock catalogue ...', end=' ')

    pos_galaxies=np.empty((N_tot,3),dtype=np.float32)
    #index: variable that go through halos (may be several galaxies in a halo)
    #i: variable that go through galaxies
    #count: find number of galaxies that lie beyond its host halo virial radius
    random.seed(seed) #this is just to check convergence on w_p(r_p)
    index=0; count=0; i=0 
    while (index<halo_mass.size):

        position=halo_pos[index]  #position of the DM halo
        radius=halo_radius[index] #radius of the DM halo

        #save the position of the central galaxy
        pos_galaxies[i]=position; i+=1

        #if halo contains satellites, save their positions
        Nsat=N_sat[index]
        if Nsat>0:
            offset=halo_offset[index] 
            length=halo_len[index]
            idss=sorted_ids[IDs[offset:offset+length]]

            #compute the radius of those particles and keep those with R<Rvir
            pos=DM_pos[idss]
            posc=pos-position

            #this is to populate correctly halos closer to box boundaries
            if np.any((position+radius>BoxSize) + (position-radius<0.0)):
                
                inside=np.where(posc[:,0]>BoxSize/2.0)[0]
                posc[inside,0]-=BoxSize
                inside=np.where(posc[:,0]<-BoxSize/2.0)[0]
                posc[inside,0]+=BoxSize

                inside=np.where(posc[:,1]>BoxSize/2.0)[0]
                posc[inside,1]-=BoxSize
                inside=np.where(posc[:,1]<-BoxSize/2.0)[0]
                posc[inside,1]+=BoxSize

                inside=np.where(posc[:,2]>BoxSize/2.0)[0]
                posc[inside,2]-=BoxSize
                inside=np.where(posc[:,2]<-BoxSize/2.0)[0]
                posc[inside,2]+=BoxSize
                
            radii=np.sqrt(posc[:,0]**2+posc[:,1]**2+posc[:,2]**2)
            inside=np.where(radii<radius)[0]
            if len(inside)<Nsat:
                problematic_cases+=1
                print('problematic case',len(inside),Nsat)
            else:
                selected=random.sample(inside,Nsat)
                pos=pos[selected]

            #aditional, not esential check. Can be comment out
            #posc=pos-position
            #if np.any((posc>BoxSize/2.0) + (posc<-BoxSize/2.0)):
            #    inside=np.where(posc[:,0]>BoxSize/2.0)[0]
            #    posc[inside,0]-=BoxSize
            #    inside=np.where(posc[:,0]<-BoxSize/2.0)[0]
            #    posc[inside,0]+=BoxSize

            #    inside=np.where(posc[:,1]>BoxSize/2.0)[0]
            #    posc[inside,1]-=BoxSize
            #    inside=np.where(posc[:,1]<-BoxSize/2.0)[0]
            #    posc[inside,1]+=BoxSize

            #    inside=np.where(posc[:,2]>BoxSize/2.0)[0]
            #    posc[inside,2]-=BoxSize
            #    inside=np.where(posc[:,2]<-BoxSize/2.0)[0]
            #    posc[inside,2]+=BoxSize
            #r_max=np.max(np.sqrt(posc[:,0]**2+posc[:,1]**2+posc[:,2]**2))
            #if r_max>radius: #check no particles beyond Rv selected 
            #    print position
            #    print radius
            #    print pos
            #    count+=1

                for j in range(Nsat):
                    pos_galaxies[i]=pos[j]; i+=1
        index+=1

    if verbose:
        print('done')
    #some final checks
    if i!=N_tot:
        print('some galaxies missing:')
        print('register',i,'galaxies out of',N_tot)
    if count>0:
        print('error:',count,'particles beyond the virial radius selected')

    return pos_galaxies
###############################################################################




##### example of use #####
"""
snapshot_fname='/data1/villa/b500p512nu0.6z99np1024tree/snapdir_017/snap_017'
groups_fname='/home/villa/data1/b500p512nu0.6z99np1024tree'
groups_number=17

### HALO CATALOGUE PARAMETERS ###
mass_criteria='t200'
min_mass=2e12 #Msun/h
max_mass=2e15 #Msun/h

### HOD PARAMETERS ###
fiducial_density=0.00111 #mean number density for galaxies with Mr<-21
M1=8e13
alpha=1.4

pos=hod(snapshot_fname,groups_fname,groups_number,min_mass,max_mass,fiducial_density,M1,alpha,mass_criteria,verbose=True)

print pos
"""
