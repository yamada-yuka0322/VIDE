import numpy as np

def cic( x, Lbox, Lboxcut = 0, Nmesh = 128, weights = None ):

 if weights == None: weights = 1
 wm = np.mean(weights)
 ws = np.mean(weights**2)

 d = np.mod(x/(Lbox+2*Lboxcut)*Nmesh,1)

 box = ([Lboxcut,Lbox+Lboxcut],[Lboxcut,Lbox+Lboxcut],[Lboxcut,Lbox+Lboxcut])

 rho = np.histogramdd(x, range = box, bins = Nmesh, weights = weights*(1-d[:,0])*(1-d[:,1])*(1-d[:,2]))[0] \
     + np.roll(np.histogramdd(x, range = box, bins = Nmesh, weights = weights*d[:,0]*(1-d[:,1])*(1-d[:,2]))[0],1,0) \
     + np.roll(np.histogramdd(x, range = box, bins = Nmesh, weights = weights*(1-d[:,0])*d[:,1]*(1-d[:,2]))[0],1,1) \
     + np.roll(np.histogramdd(x, range = box, bins = Nmesh, weights = weights*(1-d[:,0])*(1-d[:,1])*d[:,2])[0],1,2) \
     + np.roll(np.roll(np.histogramdd(x, range = box, bins = Nmesh, weights = weights*d[:,0]*d[:,1]*(1-d[:,2]))[0],1,0),1,1) \
     + np.roll(np.roll(np.histogramdd(x, range = box, bins = Nmesh, weights = weights*d[:,0]*(1-d[:,1])*d[:,2])[0],1,0),1,2) \
     + np.roll(np.roll(np.histogramdd(x, range = box, bins = Nmesh, weights = weights*(1-d[:,0])*d[:,1]*d[:,2])[0],1,1),1,2) \
     + np.roll(np.roll(np.roll(np.histogramdd(x, range = box, bins = Nmesh, weights = weights*d[:,0]*d[:,1]*d[:,2])[0],1,0),1,1),1,2)

 rho /= wm
 
 rho = rho/rho.mean() - 1.
 
 return (rho, wm, ws)
 

def powcor( d1, d2, Lbox, Nmesh = 128, Nbin = 100, scale = 'lin', cor = False ):

 # Fourier transform
 d1 = np.fft.fftn(d1)
 d2 = np.fft.fftn(d2)

 # CIC correction
 wid = np.indices(np.shape(d1)) - Nmesh/2
 #wid[np.where(wid >= Nmesh/2)] -= Nmesh
 wid = wid*np.pi/Nmesh + 1e-100
 wcic = np.prod(np.sin(wid)/wid,0)**2

 # Shell average power spectrum
 dk = 2*np.pi/Lbox
 Pk = np.conj(d1)*d2*(Lbox/Nmesh**2)**3
 Pk = np.fft.fftshift(Pk)/wcic**2

 (Nm, km, Pkm, SPkm) = shellavg(np.real(Pk), dk, Nmesh, Nbin = Nbin, xmin = 0., xmax = Nmesh*dk/2, scale = scale)
 
 # Inverse Fourier transform and shell average correlation function
 if cor == True:
  dx = Lbox/Nmesh
  Xr = np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(Pk)))*(Nmesh/Lbox)**3

  (Nmx, rm, Xrm, SXrm) = shellavg(np.real(Xr), dx, Nmesh, Nbin = Nbin, xmin = dx, xmax = 140., scale = scale)
  
  return ((Nm, km, Pkm, SPkm),(Nmx, rm, Xrm, SXrm))
 
 else: return (Nm, km, Pkm, SPkm)

 
def shellavg( f, dx, Nmesh, Nbin = 100, xmin = 0., xmax = 1., scale = 'lin' ):
 
 x = np.indices(np.shape(f)) - Nmesh/2
 #x[np.where(x >= Nmesh/2)] -= Nmesh
 x = dx*np.sqrt(np.sum(x**2,0))
 if scale == 'lin': bins = xmin+(xmax-xmin)* (np.arange(Nbin+1)/float(Nbin))
 if scale == 'log': bins = xmin*(xmax/xmin)**(np.arange(Nbin+1)/float(Nbin))

 Nm = np.histogram(x, bins = bins)[0]
 xm = np.histogram(x, bins = bins, weights = x)[0]/Nm
 fm = np.histogram(x, bins = bins, weights = f)[0]/Nm
 fs = np.sqrt((np.histogram(x, bins = bins, weights = f**2)[0]/Nm - fm**2)/(Nm-1))
 
 return (Nm, xm, fm, fs)