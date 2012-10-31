import numpy as np
import healpy as healpy
import scipy.integrate

__all__=['getSurveyProps']
   
# returns the volume and galaxy density for a given redshit slice
def getSurveyProps(maskFile, zmin, zmax, selFunMin, selFunMax, portion, selectionFuncFile=None):

  mask = healpy.read_map(maskFile)
  area = (1.*np.size(np.where(mask > 0)) / np.size(mask)) * 4.*np.pi
  zmin = zmin * 3000
  zmax = zmax * 3000
  volume = area * (zmax**3  - zmin**3) / 3

  if selectionFuncFile != None:
    selfunc = np.genfromtxt(selectionFuncFile)
    selfunc = np.array(selfunc)
    selfunc[:,0] = selfunc[:,0]/100.
    selfuncUnity = selfunc
    selfuncUnity[:,1] = 1.0
    selfuncMin = selfunc[0,0]
    selfuncMax = selfunc[-1,0]
    selfuncDx = selfunc[1,0] - selfunc[0,0]
    selfuncN = np.size(selfunc[:,0])
  
    selFunMin *= 3000
    selFunMax *= 3000
  
    selFunMin = max(selFunMin, selfuncMin)
    selFunMax = min(selFunMax, selfuncMax)
  
  
    def f(z): return selfunc[np.ceil((z-selfuncMin)/selfuncDx), 1]*z**2 
    def fTotal(z): return selfuncUnity[np.ceil((z-selfuncMin)/selfuncDx), 1]*z**2 
  
    zrange = np.linspace(selFunMin, selFunMax)
  
    nbar = scipy.integrate.quad(f, selFunMin, selFunMax)
    nbar = nbar[0]
  
    ntotal = scipy.integrate.quad(fTotal, 0.0, max(selfuncUnity[:,0]))
    #ntotal = scipy.integrate.quad(f, 0.0, max(selfunc[:,0]))
    ntotal = ntotal[0]
  
    nbar = ntotal / area / nbar

  else:
    nbar = 1.0  

  #print "PROPERTIES: ", volume, nbar
 
  return (volume, nbar)
