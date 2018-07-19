import timeit
import numpy as np
import cosmotool as ct

def myfun(N):
  f=0.001
  L = 1.0
  d=np.random.normal(size=(N,)*3)   * np.sqrt(float(N))**3 / L**3
  rho = d + f *(d*d - np.average(d*d))
  delta = (L/N)**3
  
  B = ct.bispectrum(rho * delta, 1, N, fourier=False)
  P = ct.powerspectrum(rho * delta, 1, N, fourier=False)
  PP = P[1]/P[0] / L**3 
  
  x = PP[:,None,None] * PP[None,:,None] + PP[:,None,None]*PP[None,None,:] + PP[None,:,None]*PP[None,None,:]
  
  BB = B[1]/B[0] / L**3

  y = BB/x

  np.savez("bispec_%d.npz" % N, x=x, y=y, d=d,B_nt=B[0], B_r=B[1], P_n=P[0], P=P[1], BB=BB, rho=rho, PP=PP);
  
  
#print( timeit.timeit('from __main__ import myfun; myfun(16)', number=1) )
#print( timeit.timeit('from __main__ import myfun; myfun(24)', number=1) )
print( timeit.timeit('from __main__ import myfun; myfun(32)', number=1) )
#print( timeit.timeit('from __main__ import myfun; myfun(64)', number=1) )
