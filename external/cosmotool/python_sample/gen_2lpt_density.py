import numexpr as ne
import numpy as np
import cosmotool as ct
import icgen as bic
import icgen.cosmogrowth as cg
import sys
import argparse

cosmo={'omega_M_0':0.3175, 'h':0.6711}
cosmo['omega_lambda_0']=1-cosmo['omega_M_0']
cosmo['omega_k_0'] = 0
cosmo['omega_B_0']=0.049
cosmo['SIGMA8']=0.8344
cosmo['ns']=0.9624
N0=256

doSimulation=True

astart=1.

parser=argparse.ArgumentParser(description="Generate CIC density from 2LPT")
parser.add_argument('--start', type=int, required=True)
parser.add_argument('--end', type=int, required=True)
parser.add_argument('--step', type=int, required=True)
parser.add_argument('--base', type=str, required=True)
parser.add_argument('--N', type=int, default=256)
parser.add_argument('--output', type=str, default="dcic_%d.npy")
parser.add_argument('--supersample', type=int, default=1)
parser.add_argument('--rsd', action='store_true') 
args = parser.parse_args()



for i in xrange(args.start, args.end, args.step):
  print i
#  pos,_,density,N,L,_ = bic.run_generation("/nethome/lavaux/remote/borg_2m++_128/initial_density_%d.dat" % i, 0.001, astart, cosmo, supersample=2, do_lpt2=True)
  pos,vel,density,N,L,_,_ = bic.run_generation("%s/initial_density_%d.dat" % (args.base,i), 0.001, astart, 
                                           cosmo, supersample=args.supersample, do_lpt2=True, needvel=True)

  if args.rsd:
    inv_r2 = ne.evaluate('1/sqrt(x**2+y**2+z**2)', 
        local_dict={'x':pos[0], 'y':pos[1], 'z':pos[2]})
    rsd = lambda p,v: ne.evaluate('x + (x*vx)*inv_r2 / H', 
        local_dict={'x':p, 'inv_r2':inv_r2,
         'vx':v, 'H':100.0}, out=p, casting='unsafe')

    rsd(pos[0], vel[0])
    rsd(pos[1], vel[1])
    rsd(pos[2], vel[2])

  dcic = ct.cicParticles(pos, L, args.N)
  dcic /= np.average(np.average(np.average(dcic, axis=0), axis=0), axis=0)
  dcic -= 1

  np.save(args.output % i, dcic)
