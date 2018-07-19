import numexpr as ne
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
import numpy as np
from scipy import ndimage
from scipy.special import sinc


def move_direction_new(d_theta, d_phi, theta, phi):
  cos=np.cos
  sin=np.sin
  sqrt=np.sqrt

  amplitude = sqrt(d_theta*d_theta + d_phi*d_phi);
  cos_alpha = d_theta/amplitude;
  sin_alpha = d_phi/amplitude;

  if (amplitude == 0):
    return theta,phi

  cos_d = cos(amplitude);
  sin_d = sin(amplitude);

  cos_theta = cos(theta);
  sin_theta = sin(theta);

  cos_phi = cos(phi);
  sin_phi = sin(phi);

  basis = [
     [ cos_phi * sin_theta, sin_phi * sin_theta, cos_theta ],
     [ cos_phi * cos_theta, sin_phi * cos_theta, -sin_theta ],
     [ -sin_phi, cos_phi, 0 ]
   ]

  np0 = [ cos_d, cos_alpha*sin_d, sin_alpha*sin_d ]
  np1 = [ sum([basis[j][i] * np0[j] for j in xrange(3)]) for i in xrange(3) ]
  dnp = sqrt(sum([np1[i]**2 for i in xrange(3)]))

  theta = np.arccos(np1[2]/dnp);
  phi = np.arctan2(np1[1], np1[0]) % (2*np.pi);

  return theta,phi

def move_direction_new2(delta_theta, delta_phi, theta, phi):
	cos,sin,sqrt=np.cos,np.sin,np.sqrt
	
	grad_len = sqrt(delta_theta**2 + delta_phi**2)
	if grad_len==0:
		return theta,phi

	cth0 = cos(theta)
	sth0 = sin(theta)
	
	topbottom = 1 if (theta < 0.5*np.pi) else -1
	sinc_grad_len = sinc(grad_len)

	cth = topbottom*cos(grad_len) * cth0 - sinc_grad_len*sth0*delta_theta
	sth = max(1e-10, sqrt((1.0-cth)*(1.0+cth)) )

	phi = phi + np.arcsin(delta_phi * sinc_grad_len / sth)
	theta = np.arccos(cth)
	return theta,phi

move_direction = move_direction_new

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


def build_unit_vectors(N):
  ii = np.arange(N,dtype=np.float64)/N - 0.5
  d = np.sqrt(ii[:,None,None]**2 + ii[None,:,None]**2 + ii[None,None,:]**2)
  d[N/2,N/2,N/2] = 1
  ux = ii[:,None,None] / d
  uy = ii[None,:,None] / d
  uz = ii[None,None,:] / d

  return ux,uy,uz

def compute_vcmb(l, b):
    # Motion is obtained from Tully (2007): sun_vs_LS + LS_vs_CMB
    motion = [-25.,-246.,277.];

    x = np.cos(l*np.pi/180) * np.cos(b*np.pi/180)
    y = np.sin(l*np.pi/180) * np.cos(b*np.pi/180)
    z = np.sin(b*np.pi/180)

    return x*motion[0] + y*motion[1] + z*motion[2]


def compute_vlg(l,b):

    motion = [-79,296,-36]; # [-86, 305, -33];

    x = np.cos(l*np.pi/180) * np.cos(b*np.pi/180)
    y = np.sin(l*np.pi/180) * np.cos(b*np.pi/180)
    z = np.sin(b*np.pi/180)

    return x*motion[0] + y*motion[1] + z*motion[2]


def generate_from_catalog(dmin,dmax,Nside,perturb=0.0,y=0.0,do_random=False,do_hubble=False,x=2.37,bright=-np.inf,bright_list=[],use_vlg=True,sculpt=-1):
  import progressbar as pbar
  
  cat = np.load("2m++.npy")

  cat['distance'] = cat['best_velcmb']
#  cat = cat[np.where((cat['distance']>100*dmin)*(cat['distance']<dmax*100))]

  deg2rad = np.pi/180
  Npix = 12*Nside**2
  xp,yp,zp = hp.pix2vec(Nside, np.arange(Npix))
  N2 = np.sqrt(xp**2+yp**2+zp**2)

  ksz_template = np.zeros(Npix, dtype=np.float64)
  ksz_mask = np.ones(Npix, dtype=np.uint8)
  if do_hubble:
    ksz_hubble_template = np.zeros(ksz_template.size, dtype=np.float64)

  for i in pbar.ProgressBar(maxval = cat.size, widgets=[pbar.Bar(), pbar.ETA()])(cat):
    # Skip too point sources
    if i['name'] in bright_list:
      print("Object %s is in bright list" % i['name'])
      continue

    if do_random:
      l = np.random.rand()*360
      b = np.arcsin(2*np.random.rand()-1)*180/np.pi
    else:
      l0,b0=i['gal_long'],i['gal_lat']

    l=ne.evaluate('l0*deg2rad')
    b=ne.evaluate('b0*deg2rad')

    dtheta,dphi = np.random.randn(2)*perturb
    theta,l=move_direction(dtheta,dphi,0.5*np.pi - b, l)
	
    b = 0.5*np.pi-theta

    x0 = np.cos(l)*np.cos(b)
    y0 = np.sin(l)*np.cos(b)
    z0 = np.sin(b)

    if use_vlg:
      vlg = i['best_velcmb'] - compute_vcmb(l0, b0) + compute_vlg(l0, b0)
      DA = vlg/100
    else:
      DA = i['best_velcmb'] / 100

    if DA < dmin or DA > dmax:
      continue

    Lgal = DA**2*10**(0.4*(tmpp_cat['Msun']-i['K2MRS']+25))

    M_K=i['K2MRS']-5*np.log10(DA)-25
    # Skip too bright galaxies
    if M_K < bright:
      continue

    profiler = ksz.KSZ_Isothermal(Lgal, x, y=y, sculpt=sculpt)

    idx0 = hp.query_disc(Nside, (x0,y0,z0), 3*profiler.rGalaxy/DA)

    xp1 = xp[idx0]
    yp1 = yp[idx0]
    zp1 = zp[idx0]
    N2_1 = N2[idx0]

    cos_theta = ne.evaluate('(x0*xp1+y0*yp1+z0*zp1)/(sqrt(x0**2+y0**2+z0**2)*(N2_1))')

    idx,idx_masked,m = profiler.projected_profile(cos_theta, DA)
    idx = idx0[idx]
    idx_masked = idx0[idx_masked]
    ksz_template[idx] += m
    ksz_mask[idx_masked] = 0
    if do_hubble:
      ksz_hubble_template[idx] += m*DA

  ne.evaluate('ksz_template*ksz_normalization', out=ksz_template)

  result =ksz_template, ksz_mask
  if do_hubble:
    ne.evaluate('ksz_hubble_template*ksz_normalization', out=ksz_hubble_template)
    return result + ( ksz_hubble_template,)
  else:
    return result

def get_args():
  parser=argparse.ArgumentParser(description="Generate Skymaps from CIC maps")
  parser.add_argument('--Nside', type=int, default=128)
  parser.add_argument('--minval', type=float, default=-0.5)
  parser.add_argument('--maxval', type=float, default=0.5)
  parser.add_argument('--depth_min', type=float, default=10)
  parser.add_argument('--depth_max', type=float, default=60)
  parser.add_argument('--ksz_map', type=str, required=True)
  parser.add_argument('--base_fig', type=str, default="kszfig.png")
  parser.add_argument('--build_dipole', action='store_true')
  parser.add_argument('--degrade', type=int, default=-1)
  parser.add_argument('--y',type=float,default=0.0)
  parser.add_argument('--x',type=float,default=2.37)
  parser.add_argument('--random', action='store_true')
  parser.add_argument('--perturb', type=float, default=0)
  parser.add_argument('--hubble_monopole', action='store_true')
  parser.add_argument('--remove_bright', type=float, default=-np.inf)
  parser.add_argument('--bright_file', type=str)
  parser.add_argument('--lg', action='store_true')
  parser.add_argument('--sculpt_beam', type=float, default=-1)
  return parser.parse_args()

def main():

  args = get_args()

  ff=plt.figure(1)
  plt.clf()
  v=[]

  print("Generating map...")

  with open("crap.txt", mode="r") as f:
    bright_list = [l.split('#')[0].strip(" \t\n\r") for l in f]

  if args.bright_file:
    with open(args.bright_file, mode="r") as f:
      idx_name = f.readline().split(',').index('name_2')
      bright_list = bright_list + [l.split(',')[idx_name] for l in f]

  print("Built bright point source list: " + repr(bright_list))
  r = generate_from_catalog(args.depth_min,args.depth_max,args.Nside,perturb=args.perturb,y=args.y,do_random=args.random,do_hubble=args.hubble_monopole,x=args.x,bright=args.remove_bright,use_vlg=args.lg,bright_list=bright_list,sculpt=args.sculpt_beam)
  hubble_map = None
  if args.hubble_monopole:
    proj,mask,hubble_map = r
  else:
    proj,mask = r

  if args.degrade > 0:
    proj *= mask
    proj = hp.ud_grade(proj, nside_out=args.degrade)
    if hubble_map is not None:
      hubble_map *= mask
      hubble_map = hp.ud_grade(hubble_map, nside_out=args.degrade)
    mask = hp.ud_grade(mask, nside_out=args.degrade)
    Nside = args.degrade
  else:
    Nside = args.Nside

  hp.write_map(args.ksz_map + ".fits", proj)
  hp.write_map(args.ksz_map + "_mask.fits", mask)

  if args.build_dipole:
    x,y,z=hp.pix2vec(Nside, np.arange(hp.nside2npix(Nside)))
    hp.write_map(args.ksz_map + "_x.fits", proj*x)
    hp.write_map(args.ksz_map + "_y.fits", proj*y)
    hp.write_map(args.ksz_map + "_z.fits", proj*z)

  if args.hubble_monopole:
    hp.write_map(args.ksz_map + "_hubble.fits", hubble_map)

  hp.mollview(proj*100*1e6, fig=1, coord='GG', cmap=plt.cm.coolwarm, title='', min=args.minval, 
              max=args.maxval)

  ff.savefig(args.base_fig)
  
main()
