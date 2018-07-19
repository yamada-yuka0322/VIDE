import healpy as hp
import numpy as np
import cosmotool as ct
import h5py as h5
from matplotlib import pyplot as plt

L=600.
Nside=128

INDATA="/nethome/lavaux/Copy/PlusSimulation/BORG/Input_Data/2m++.npy"
tmpp = np.load(INDATA)

def build_sky_proj(density, dmax=60.,dmin=0):

  N = density.shape[0]
  ix = (np.arange(N)-0.5)*L/N - 0.5 * L


  dist2 = (ix[:,None,None]**2 + ix[None,:,None]**2 + ix[None,None,:]**2)

  flux = density.transpose().astype(ct.DTYPE) # / dist2
  dmax=N*dmax/L
  dmin=N*dmin/L
  projsky1 = ct.spherical_projection(Nside, flux, dmin, dmax, integrator_id=1)
#  projsky0 = ct.spherical_projection(Nside, flux, 0, 52, integrator_id=0)

  return projsky1*L/N#,projsky0

l,b = tmpp['gal_long'],tmpp['gal_lat']

l = np.radians(l)
b = np.pi/2 - np.radians(b)

dcmb = tmpp['velcmb']/100.

idx = np.where((dcmb>10)*(dcmb<60))

plt.figure(1)
plt.clf()
if True:
  with h5.File("fields.h5", mode="r") as f:
    d = f["density"][:].transpose()
    d /= np.average(np.average(np.average(d,axis=0),axis=0),axis=0)
    proj = build_sky_proj(d, dmin=10,dmax=60.)
  proj0 = proj1 = proj
else:
  d = np.load("icgen/dcic0.npy")
  proj0 = build_sky_proj(1+d, dmin=10,dmax=60.)
  d = np.load("icgen/dcic1.npy")
  proj1 = build_sky_proj(1+d, dmin=10,dmax=60.)

hp.mollview(proj0, fig=1, coord='CG', max=60, cmap=plt.cm.coolwarm)
hp.projscatter(b[idx], l[idx], lw=0, color='g', s=5.0, alpha=0.8)

plt.figure(2)
plt.clf()
hp.mollview(proj1, fig=2, coord='CG', max=60, cmap=plt.cm.coolwarm)
hp.projscatter(b[idx], l[idx], lw=0, color='g', s=5.0, alpha=0.8)
