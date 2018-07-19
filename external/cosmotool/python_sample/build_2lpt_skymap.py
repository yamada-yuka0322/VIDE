import matplotlib
matplotlib.use('Agg')
import healpy as hp
import numpy as np
import cosmotool as ct
import argparse
import h5py as h5
from matplotlib import pyplot as plt

parser=argparse.ArgumentParser(description="Generate Skymaps from CIC maps")
parser.add_argument('--boxsize', type=float, required=True)
parser.add_argument('--Nside', type=int, default=128)
parser.add_argument('--base_cic', type=str, required=True)
parser.add_argument('--base_fig', type=str, required=True)
parser.add_argument('--start', type=int, required=True)
parser.add_argument('--end', type=int, required=True)
parser.add_argument('--step', type=int, required=True)
parser.add_argument('--minval', type=float, default=0)
parser.add_argument('--maxval', type=float, default=4)
parser.add_argument('--depth_min', type=float, default=10)
parser.add_argument('--depth_max', type=float, default=60)
parser.add_argument('--iid', type=int, default=0)
parser.add_argument('--proj_cat', type=bool, default=False)
args = parser.parse_args()

#INDATA="/nethome/lavaux/Copy/PlusSimulation/BORG/Input_Data/2m++.npy"
INDATA="2m++.npy"
tmpp = np.load(INDATA)

L = args.boxsize
Nside = args.Nside

def build_sky_proj(density, dmax=60.,dmin=0,iid=0):

  N = density.shape[0]
  ix = (np.arange(N)-0.5)*L/N - 0.5 * L

#  dist2 = (ix[:,None,None]**2 + ix[None,:,None]**2 + ix[None,None,:]**2)

  flux = density.transpose().astype(ct.DTYPE) # / dist2
  dmax=N*dmax/L
  dmin=N*dmin/L
  if iid == 0:
    shifter = np.array([0.5,0.5,0.5])
  else:
    shifter = np.array([0.,0.,0.])
    
  projsky1 = ct.spherical_projection(Nside, flux, dmin, dmax, integrator_id=iid, shifter=shifter)

  return projsky1*L/N

l,b = tmpp['gal_long'],tmpp['gal_lat']

l = np.radians(l)
b = np.pi/2 - np.radians(b)

dcmb = tmpp['velcmb']/100.

idx = np.where((dcmb>args.depth_min)*(dcmb<args.depth_max))


for i in xrange(args.start,args.end,args.step):
  ff=plt.figure(1)
  plt.clf()
  d = np.load(args.base_cic % i)
  proj = build_sky_proj(1+d, dmin=args.depth_min,dmax=args.depth_max,iid=args.iid)
  proj /= (args.depth_max-args.depth_min)

  hp.write_map("skymaps/proj_map_%d.fits" % i, proj)

  print proj.min(), proj.max()
  hp.mollview(proj, fig=1, coord='CG', cmap=plt.cm.copper, title='Sample %d' % i, min=args.minval, max=args.maxval)
  if args.proj_cat:
    hp.projscatter(b[idx], l[idx], lw=0, color=[0.1,0.8,0.8], s=2.0, alpha=0.7)

  ff.savefig(args.base_fig % i)
