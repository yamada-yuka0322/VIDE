import healpy as hp
import numpy as np
import cosmotool as ct
import argparse
import h5py as h5
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import ksz
from ksz.constants import *
from cosmotool import interp3d

def wrapper_impulsion(f):

  class _Wrapper(object):
    def __init__(self):
      pass
      
    def __getitem__(self,direction):
    
      if 'velocity' in f:
        return f['velocity'][:,:,:,direction]
      
      n = "p%d" % direction
      return f[n]

  return _Wrapper()

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

def build_unit_vectors(N):
  ii = np.arange(N)*L/N - 0.5*L
  d = np.sqrt(ii[:,None,None]**2 + ii[None,:,None]**2 + ii[None,None,:]**2)
  d[N/2,N/2,N/2] = 1
  ux = ii[:,None,None] / d
  uy = ii[None,:,None] / d
  uz = ii[None,None,:] / d

  return ux,uy,uz

def build_radial_v(v):
  N = v[0].shape[0]
  u = build_unit_vectors(N)
  vr = v[0] * u[2]
  vr += v[1] * u[1]
  vr += v[2] * u[0]

  return vr.transpose()

def generate_from_catalog(vfield,Boxsize,dmin,dmax):
  import progressbar as pbar
  
  cat = np.load("2m++.npy")


  cat['distance'] = cat['best_velcmb']
  cat = cat[np.where((cat['distance']>100*dmin)*(cat['distance']<dmax*100))]

  deg2rad = np.pi/180
  Npix = 12*Nside**2
  xp,yp,zp = hp.pix2vec(Nside, np.arange(Npix))
  N2 = np.sqrt(xp**2+yp**2+zp**2)

  ksz_template = np.zeros(Npix, dtype=np.float64)
  ksz_mask = np.zeros(Npix, dtype=np.uint8)

  pb = pbar.ProgressBar(maxval = cat.size, widgets=[pbar.Bar(), pbar.ETA()]).start()

  for k,i in np.ndenumerate(cat):
    pb.update(k[0])
    l,b=i['gal_long'],i['gal_lat']
    ra,dec=i['ra'],i['dec']

    l *= deg2rad
    b *= deg2rad

    ra *= deg2rad
    dec *= deg2rad

#    x0 = np.cos(l)*np.cos(b)
#    y0 = np.sin(l)*np.cos(b)
#    z0 = np.sin(b)

    x0 =  xra = np.cos(ra)*np.cos(dec)
    y0 = yra = np.sin(ra)*np.cos(dec)
    z0 = zra = np.sin(dec)
    
    DA =i['distance']/100
    Lgal = DA**2*10**(0.4*(tmpp_cat['Msun']-i['K2MRS']+25))

    profiler = ksz.KSZ_Isothermal(Lgal, 2.37)

    idx0 = hp.query_disc(Nside, (x0,y0,z0), 3*profiler.rGalaxy/DA)

    vr = interp3d(DA * xra, DA * yra, DA * zra, vfield, Boxsize)

    xp1 = xp[idx0]
    yp1 = yp[idx0]
    zp1 = zp[idx0]
    N2_1 = N2[idx0]

    cos_theta = x0*xp1+y0*yp1+z0*zp1
    cos_theta /= np.sqrt(x0**2+y0**2+z0**2)*(N2_1)

    idx,idx_masked,m = profiler.projected_profile(cos_theta, DA)
    idx = idx0[idx]
    idx_masked = idx0[idx_masked]
    ksz_template[idx] += vr * m
    ksz_mask[idx_masked] = 0

  pb.finish()

  return ksz_template, ksz_mask

for i in xrange(args.start,args.end,args.step):
  ff=plt.figure(1)
  plt.clf()
  v=[]
  fname = args.base_h5 % i
  if False:
    print("Opening %s..." % fname)
    with h5.File(fname, mode="r") as f:
      p = wrapper_impulsion(f)
      for j in xrange(3):
        v.append(p[j] / f['density'][:])

    print("Building radial velocities...")

  # Rescale by Omega_b / Omega_m
    vr = build_radial_v(v)
#  _,_,vr = build_unit_vectors(128)
#  vr *= 1000 * 500
    vr *= ksz_normalization*1e6
    del v
#    np.save("vr.npy", vr)
  else:
    print("Loading vrs...")
    vr = np.load("vr.npy")

  print("Generating map...")

  proj,mask = generate_from_catalog(vr,args.boxsize,args.depth_min,args.depth_max)

  hp.write_map(args.ksz_map % i, proj)
  hp.write_map((args.ksz_map % i) + "_mask", mask)

  hp.mollview(proj, fig=1, coord='CG', cmap=plt.cm.coolwarm, title='Sample %d' % i, min=args.minval, 
              max=args.maxval)

  ff.savefig(args.base_fig % i)
