import numpy as np
import cosmotool as ct
import borgicgen as bic
import cosmogrowth as cg
import sys

cosmo={'omega_M_0':0.3175, 'h':0.6711}
cosmo['omega_lambda_0']=1-cosmo['omega_M_0']
cosmo['omega_k_0'] = 0
cosmo['omega_B_0']=0.049
cosmo['SIGMA8']=0.8344
cosmo['ns']=0.9624
N0=256

doSimulation=False
simShift=False

snap_id=int(sys.argv[1])
astart=1/100.

if doSimulation:
  s = ct.loadRamsesAll("/nethome/lavaux/remote2/borgsim3/", snap_id, doublePrecision=True)
  astart=s.getTime()
  L = s.getBoxsize()

  p = s.getPositions()
  Nsim = int( np.round( p[0].size**(1./3)) )
  print("Nsim = %d" % Nsim)

  if simShift:
    p = [(q-0.5*L/Nsim)%L for q in p]

  dsim = ct.cicParticles(p[::-1], L, N0)
  dsim /= np.average(np.average(np.average(dsim, axis=0), axis=0), axis=0)
  dsim -= 1

  dsim_hat = np.fft.rfftn(dsim)*(L/N0)**3
  Psim, bsim = bic.bin_power(np.abs(dsim_hat)**2/L**3, L, range=(0,1.), bins=150)

pos,_,density,N,L,_,_ = bic.run_generation("initial_density_1872.dat", 0.001, astart, cosmo, supersample=1, do_lpt2=False, supergenerate=2)

dcic = ct.cicParticles(pos, L, N0)
dcic /= np.average(np.average(np.average(dcic, axis=0), axis=0), axis=0)
dcic -= 1

dcic_hat = np.fft.rfftn(dcic)*(L/N0)**3
dens_hat = np.fft.rfftn(density)*(L/N0)**3

Pcic, bcic = bic.bin_power(np.abs(dcic_hat)**2/L**3, L, range=(0,4.), bins=150)
Pdens, bdens = bic.bin_power(np.abs(dens_hat)**2/L**3, L, range=(0,4.), bins=150)

cgrowth = cg.CosmoGrowth(**cosmo)
D1 = cgrowth.D(astart)
D1_0 = D1/cgrowth.D(1)#0.001)

Pref, bref = bic.compute_ref_power(L, N0, cosmo, range=(0,4.), bins=150)

Pcic /= D1_0**2

#borg_evolved = ct.read_borg_vol("final_density_1380.dat")
#dborg_hat = np.fft.rfftn(borg_evolved.density)*L**3/borg_evolved.density.size

#Pborg, bborg = bic.bin_power(np.abs(dborg_hat)**2/L**3, L, range=(0,1.),bins=150)
