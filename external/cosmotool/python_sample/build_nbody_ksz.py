import healpy as hp
import numpy as np
import cosmotool as ct
import argparse
import h5py as h5
from matplotlib import pyplot as plt

Mpc=3.08e22
rhoc = 1.8864883524081933e-26  # m^(-3)
sigmaT = 6.6524e-29
mp = 1.6726e-27
lightspeed = 299792458.
v_unit = 1e3  # Unit of 1 km/s
T_cmb=2.725
h = 0.71
Y = 0.245  #The Helium abundance
Omega_matter = 0.26
Omega_baryon=0.0445

G=6.67e-11
MassSun=2e30
frac_electron = 1.0 # Hmmmm
frac_gas_galaxy = 0.14
mu = 1/(1-0.5*Y)

baryon_fraction = Omega_baryon / Omega_matter

ksz_normalization = T_cmb*sigmaT*v_unit/(lightspeed*mu*mp) * baryon_fraction
rho_mean_matter = Omega_matter * (3*(100e3/Mpc)**2/(8*np.pi*G)) 

parser=argparse.ArgumentParser(description="Generate Skymaps from CIC maps")
parser.add_argument('--boxsize', type=float, required=True)
parser.add_argument('--Nside', type=int, default=128)
parser.add_argument('--base_h5', type=str, required=True)
parser.add_argument('--base_fig', type=str, required=True)
parser.add_argument('--start', type=int, required=True)
parser.add_argument('--end', type=int, required=True)
parser.add_argument('--step', type=int, required=True)
parser.add_argument('--minval', type=float, default=-0.5)
parser.add_argument('--maxval', type=float, default=0.5)
parser.add_argument('--depth_min', type=float, default=10)
parser.add_argument('--depth_max', type=float, default=60)
parser.add_argument('--iid', type=int, default=0)
parser.add_argument('--ksz_map', type=str, required=True)
args = parser.parse_args()

L = args.boxsize
Nside = args.Nside

def build_sky_proj(density, dmax=60.,dmin=0,iid=0):

  N = density.shape[0]
  ix = (np.arange(N)-0.5)*L/N - 0.5 * L


  dist2 = (ix[:,None,None]**2 + ix[None,:,None]**2 + ix[None,None,:]**2)

  flux = density.transpose().astype(ct.DTYPE) # / dist2
  dmax=N*dmax/L
  dmin=N*dmin/L
  if iid == 0:
    shifter = np.array([0.5,0.5,0.5])
  else:
    shifter = np.array([0.,0.,0.])
    
  projsky1 = ct.spherical_projection(Nside, flux, dmin, dmax, integrator_id=iid, shifter=shifter)

  return projsky1*L/N

def build_unit_vectors(N):
  ii = np.arange(N)*L/N - 0.5*L
  d = np.sqrt(ii[:,None,None]**2 + ii[None,:,None]**2 + ii[None,None,:]**2)
  d[N/2,N/2,N/2] = 1
  ux = ii[:,None,None] / d
  uy = ii[None,:,None] / d
  uz = ii[None,None,:] / d

  return ux,uy,uz

for i in xrange(args.start,args.end,args.step):
  ff=plt.figure(1)
  plt.clf()
  with h5.File(args.base_h5 % i, mode="r") as f:
    p = f['velocity'][:]
    davg = np.average(np.average(np.average(f['density'][:],axis=0),axis=0),axis=0)
    p /= davg # Now we have momentum scaled to the mean density

  # Rescale by Omega_b / Omega_m
  p = p.astype(np.float64)
  print p.max(), p.min(), ksz_normalization, rho_mean_matter
  p *= -ksz_normalization*rho_mean_matter*1e6
  
  ux,uy,uz = build_unit_vectors(p.shape[0])
  pr = p[:,:,:,0] * ux + p[:,:,:,1] * uy + p[:,:,:,2] * uz 
  print p.max(), p.min()
  print pr.max()*Mpc, pr.min()*Mpc
  
  @ct.timeit_quiet
  def run_proj():
    return build_sky_proj(pr*Mpc, dmin=args.depth_min,dmax=args.depth_max,iid=args.iid)

  run_proj()

  hp.write_map(args.ksz_map % i, proj)

  hp.mollview(proj, fig=1, coord='CG', cmap=plt.cm.coolwarm, title='Sample %d' % i, min=args.minval, 
              max=args.maxval)

  ff.savefig(args.base_fig % i)
