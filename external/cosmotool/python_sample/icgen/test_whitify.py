import numpy as np
import cosmotool as ct
import borgicgen as bic
from matplotlib import pyplot as plt

cosmo={'omega_M_0':0.3175, 'h':0.6711}
cosmo['omega_lambda_0']=1-cosmo['omega_M_0']
cosmo['omega_k_0'] = 0
cosmo['omega_B_0']=0.049
cosmo['SIGMA8']=0.8344
cosmo['ns']=0.9624

zstart=50
astart=1/(1.+zstart)
halfPixelShift=False

posx,vel,density,N,L,a_ic,cosmo = bic.run_generation("initial_condition_borg.dat", 0.001, astart, cosmo, supersample=1, shiftPixel=halfPixelShift, do_lpt2=False)

w1 = bic.whitify(density, L, cosmo, supergenerate=1)
w2 = bic.whitify(density, L, cosmo, supergenerate=2)

N = w1.shape[0]
Ns = w2.shape[0]

w1_hat = np.fft.rfftn(w1)*(L/N)**3
w2_hat = np.fft.rfftn(w2)*(L/Ns)**3

P1, b1, dev1 = bic.bin_power(np.abs(w1_hat)**2, L, range=(0,3),bins=150,dev=True)
P2, b2, dev2 = bic.bin_power(np.abs(w2_hat)**2, L, range=(0,3),bins=150,dev=True)

fig = plt.figure(1)
fig.clf()
plt.fill_between(b1, P1*(1-dev1), P1*(1+dev1), label='Supergen=1', color='b')
plt.fill_between(b2, P2*(1-dev2), P2*(1+dev2), label='Supergen=2', color='g', alpha=0.5)
ax = plt.gca()
ax.set_xscale('log')
plt.ylim(0.5,1.5)
plt.xlim(1e-2,4)
plt.axhline(1.0, color='red', lw=4.0)
plt.legend()
plt.show()
