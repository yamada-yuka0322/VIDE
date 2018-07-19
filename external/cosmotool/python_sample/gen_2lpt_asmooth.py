import os
import h5py as h5
import numpy as np
import cosmotool as ct
import icgen as bic
import icgen.cosmogrowth as cg
import sys
import argparse

#ADAPT_SMOOTH="/home/bergeron1NS/lavaux/Software/cosmotool/build/sample/simple3DFilter"
ADAPT_SMOOTH="/home/guilhem/PROJECTS/cosmotool/build/sample/simple3DFilter"
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
parser.add_argument('--output', type=str, default="fields_%d.h5")
parser.add_argument('--supersample', type=int, default=1)
args = parser.parse_args()



for i in [4629]:#xrange(args.start, args.end, args.step):
  print i

  pos,vel,density,N,L,_,_ = bic.run_generation("%s/initial_density_%d.dat" % (args.base,i), 0.001, astart, 
                                           cosmo, supersample=args.supersample, do_lpt2=True)

  q = pos + vel + [np.ones(vel[0].shape[0])]

  with h5.File("particles.h5", mode="w") as f:
    f.create_dataset("particles", data=np.array(q).transpose())
  
  os.system(ADAPT_SMOOTH + " %s %lg %lg %d %lg %lg %lg" % ("particles.h5", 3000000, L, args.N, 0, 0, 0))
  os.rename("fields.h5", args.output % i)
