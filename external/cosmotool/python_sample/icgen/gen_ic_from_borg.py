import pyfftw
import numpy as np
import cosmotool as ct
import borgicgen as bic
import pickle

with file("wisdom") as f:
  pyfftw.import_wisdom(pickle.load(f))

cosmo={'omega_M_0':0.3175, 'h':0.6711}
cosmo['omega_lambda_0']=1-cosmo['omega_M_0']
cosmo['omega_k_0'] = 0
cosmo['omega_B_0']=0.049
cosmo['SIGMA8']=0.8344
cosmo['ns']=0.9624

supergen=1
zstart=99
astart=1/(1.+zstart)
halfPixelShift=False
zero_fill=False

if __name__=="__main__":
    bic.write_icfiles(*bic.run_generation("initial_density_1872.dat", 0.001, astart, cosmo, supersample=1, shiftPixel=halfPixelShift, do_lpt2=False, supergenerate=supergen), supergenerate=1, zero_fill=zero_fill)
