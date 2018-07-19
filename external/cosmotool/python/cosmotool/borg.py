### 
### BORG code is from J. Jasche
###
import io
import numpy as np
from numpy import *
import os.path
import array
import glob

class BorgVolume(object):

    def __init__(self, density, ranges):
        self.density = density
        self.ranges = ranges

def build_filelist(fdir):
    #builds list of all borg density fields which may be distributed over several directories
    
    fname_0=glob.glob(fdir[0]+'initial_density_*')
    fname_1=glob.glob(fdir[0]+'final_density_*')
       
    fdir=fdir[1:] #eliminate first element
       
    for fd in fdir:
      fname_0=fname_0+glob.glob(fd+'initial_density_*')
      fname_1=fname_1+glob.glob(fd+'final_density_*')
	   
    return fname_0, fname_1

def read_borg_vol(BORGFILE):
    """ Reading routine for BORG data
    """

    openfile=open(BORGFILE,'rb')

    period=0
    N0=0
    N1=0
    N2=0

    xmin=0
    xmax=0

    ymin=0
    ymax=0

    zmin=0
    zmax=0
    nlines=0

    while True:
      line=openfile.readline()
      s=line.rstrip('\n')
      r=s.rsplit(' ')
    
      if size(r)==5 : 
        if r[0] =="define":
          if r[1]=="Lattice" :
            N0=int(r[2])
            N1=int(r[3])
            N2=int(r[4])
            
          if size(r)==11 :
            if r[4] =="BoundingBox":
              xmin=float(r[5])
              xmax=float(r[6])
              ymin=float(r[7])
              ymax=float(r[8])
              zmin=float(r[9])
              zmax=float(r[10].rstrip(','))
              
          if r[0]=='@1': break
  
    ranges=[]
    ranges.append(xmin)
    ranges.append(xmax)
    ranges.append(ymin)
    ranges.append(ymax)
    ranges.append(zmin)
    ranges.append(zmax)
    
    #now read data
    data=np.fromfile(openfile, '>f4')
    data=data.reshape(N2,N0,N1)
    return BorgVolume(data,ranges)

def read_spec( fname ):
    """ Reading routine for ARES spectrum samples
    """
    x,y=np.loadtxt( fname ,usecols=(0,1),unpack=True)

    return x , y
    
        
def read_bias_nmean( fname ):
    """ Reading routine for ARES bias data
    """
    x,b0,b1,nmean=np.loadtxt( fname ,usecols=(0,1,2,3),unpack=True)

    return x , b0, b1, nmean

def read_nmean( fname ):
    """ Reading routine for BORG bias data
    """
    x,nmean=np.loadtxt( fname ,usecols=(0,1),unpack=True)

    return x, nmean    
    
def get_grid_values(xx,data, ranges):
    """ return values at grid positions
    """
    xmin=ranges[0]
    xmax=ranges[1]
    
    ymin=ranges[2]
    ymax=ranges[3]
    
    zmin=ranges[4]
    zmax=ranges[5]
    
    Lx= xmax-xmin
    Ly= ymax-ymin
    Lz= zmax-zmin
        
    Nx=shape(data)[0]
    Ny=shape(data)[1]
    Nz=shape(data)[2]
    
    dx=Lx/float(Nx)
    dy=Ly/float(Ny)
    dz=Lz/float(Nz)
    
    idx=(xx[:,0]-xmin)/dx
    idy=(xx[:,1]-ymin)/dz
    idz=(xx[:,2]-zmin)/dy
    
    idx=idx.astype(int)
    idy=idy.astype(int)
    idz=idz.astype(int)
    
    idflag=np.where( (idx>-1)*(idx<Nx)*(idy>-1)*(idy<Ny)*(idz>-1)*(idz<Nz) )
        
    flag=[False]*len(xx)
    vals=[-999.]*len(xx)
    
    flag=np.array(flag)
    vals=np.array(vals)
    flag[idflag]=True
    vals[idflag]=data[idx[idflag],idy[idflag],idz[idflag]]
        
    return vals,flag    
    
def get_mean_density(fdir, smin, step):
    """ estimate ensemble mean
    """
    import progressbar as pb

    print('-'*60)
    print('Get 3D ensemble mean density field')
    print('-'*60)

    fname0 = fdir + 'initial_density_'+str(0)+'.dat'
    fname1 = fdir + 'final_density_'+str(0)+'.dat'
        
    MEAN0,ranges=read_borg_vol(fname0)
    MEAN0=MEAN0*0.;
    VAR0=copy(MEAN0)
    MEAN1=copy(MEAN0)
    VAR1=copy(MEAN0)
    norm=0.
    
    idat=smin
    
    fname0 = fdir + 'initial_density_'+str(idat)+'.dat'
    fname1 = fdir + 'final_density_'+str(idat)+'.dat'
        #and (idat<smin+1000)
    while((os.path.exists(fname0))):
      auxdata0,auxranges0=read_borg_vol(fname0)
      auxdata1,auxranges1=read_borg_vol(fname1)
      auxx0=auxdata0
      auxx1=auxdata1
      MEAN0+=auxx0
      VAR0+=auxx0**2
      MEAN1+=auxx1
      VAR1+=auxx1**2

      norm+=1
      idat+=step
      fname0 = fdir + 'initial_density_'+str(idat)+'.dat'
      fname1 = fdir + 'final_density_'+str(idat)+'.dat'
      del auxranges0
      del auxdata0
      del auxranges1
      del auxdata1
	    
    MEAN0/=norm
    VAR0/=norm
    VAR0-=MEAN0**2
    VAR0=sqrt(fabs(VAR0))
    
    MEAN1/=norm
    VAR1/=norm
    VAR1-=MEAN1**2
    VAR1=sqrt(fabs(VAR1))
    
    
    return MEAN0,VAR0,MEAN1,VAR1,ranges

def get_mean_density_fdir(fdir,init,steps):
    """ estimate ensemble mean
    """
    import progressbar as pb

    print('-'*60)
    print('Get 3D ensemble mean density field')
    print('-'*60)

    fname0,fname1=build_filelist(fdir)
    
    fname0=fname0[init::steps]
    fname1=fname1[init::steps]
    
    borg=read_borg_vol(fname0[0])
    MEAN0 = borg.density
    RANGES0 = borg.ranges
    MEAN0=MEAN0*0.;
    VAR0=copy(MEAN0)
    MEAN1=copy(MEAN0)
    VAR1=copy(MEAN0)
    norm0=0.
    norm1=0.
    
    for fn in pb.ProgressBar(len(fname0))(fname0):
      auxborg=read_borg_vol(fn)
      auxdata0 = auxborg.density
      MEAN0+=auxdata0
      VAR0+=auxdata0**2.
      norm0+=1.
      del auxdata0
      del auxborg

    for fn in pb.ProgressBar(len(fname1))(fname1):
      auxborg1=read_borg_vol(fn)
      auxdata1 = auxborg1.density
      MEAN1+=auxdata1
      VAR1+=auxdata1**2.
      norm1+=1.
      del auxdata1
      del auxborg1
	
    MEAN0/=norm0
    VAR0/=norm0
    VAR0-=MEAN0**2
    VAR0=sqrt(fabs(VAR0))
    
    MEAN1/=norm1
    VAR1/=norm1
    VAR1-=MEAN1**2
    VAR1=sqrt(fabs(VAR1))
    
    
    return MEAN0,VAR0,MEAN1,VAR1,ranges    
    
    
