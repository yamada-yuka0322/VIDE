import cosmotool as ct
import numpy as np
import cosmolopy as cpy
from cosmogrowth import *
import borgadaptor as ba

@ct.timeit
def gen_posgrid(N, L, delta=1, dtype=np.float32):
    """ Generate an ordered lagrangian grid"""

    ix = (np.arange(N)*(L/N*delta)).astype(dtype)

    x = ix[:,None,None].repeat(N, axis=1).repeat(N, axis=2)
    y = ix[None,:,None].repeat(N, axis=0).repeat(N, axis=2)
    z = ix[None,None,:].repeat(N, axis=0).repeat(N, axis=1)

    return x.reshape((x.size,)), y.reshape((y.size,)), z.reshape((z.size,))

def bin_power(P, L, bins=20, range=(0,1.), dev=False):

  N = P.shape[0]
  ik = np.fft.fftfreq(N, d=L/N)*2*np.pi

  k = np.sqrt(ik[:,None,None]**2 + ik[None,:,None]**2 + ik[None,None,:(N/2+1)]**2)

  H,b = np.histogram(k, bins=bins, range=range)
  Hw,b = np.histogram(k, bins=bins, weights=P, range=range)
  
  if dev:
    return Hw/(H-1), 0.5*(b[1:]+b[0:bins]), 1.0/np.sqrt(H)
  else:
    return Hw/(H-1), 0.5*(b[1:]+b[0:bins])

def compute_power_from_borg(input_borg, a_borg, cosmo, bins=10, range=(0,1)):
  borg_vol = ct.read_borg_vol(input_borg)
  N = borg_vol.density.shape[0]

  cgrowth = CosmoGrowth(**cosmo)
  D1 = cgrowth.D(1)
  D1_0 = D1/cgrowth.D(a_borg)
  print("D1_0=%lg" % D1_0)

  density_hat, L = ba.half_pixel_shift(borg_vol)

  return bin_power(D1_0**2*np.abs(density_hat)**2/L**3, L, bins=bins, range=range)
  
def compute_ref_power(L, N, cosmo, bins=10, range=(0,1), func='HU_WIGGLES'):
  ik = np.fft.fftfreq(N, d=L/N)*2*np.pi

  k = np.sqrt(ik[:,None,None]**2 + ik[None,:,None]**2 + ik[None,None,:(N/2+1)]**2)
  p = ct.CosmologyPower(**cosmo)
  p.setFunction(func)
  p.normalize(cosmo['SIGMA8'])

  return bin_power(p.compute(k)*cosmo['h']**3, L, bins=bins, range=range)


def do_supergenerate(density, density_out=None, mulfac=None,zero_fill=False,Pk=None,L=None,h=None):
 
  N = density.shape[0]
  
  if density_out is None:
    assert mulfac is not None
    Ns = mulfac*N
    density_out = np.zeros((Ns,Ns,Ns/2+1), dtype=np.complex128)
    density_out[:] = np.nan
  elif mulfac is None:
    mulfac = density_out.shape[0] / N
    Ns = density_out.shape[0]
    assert (density_out.shape[0] % N) == 0
  	
  assert len(density_out.shape) == 3
  assert density_out.shape[0] == density_out.shape[1]
  assert density_out.shape[2] == (density_out.shape[0]/2+1)
  
  hN = N/2
  density_out[:hN,               :hN, :hN+1] = density[:hN, :hN, :]
  density_out[:hN,        (Ns-hN):Ns, :hN+1] = density[:hN, hN:, :]
  density_out[(Ns-hN):Ns, (Ns-hN):Ns, :hN+1] = density[hN:, hN:, :]
  density_out[(Ns-hN):Ns,        :hN, :hN+1] = density[hN:, :hN, :]

  if mulfac > 1:

    cond=np.isnan(density_out)
    if zero_fill:
      density_out[cond] = 0
    else:

      if Pk is not None:
        assert L is not None and h is not None

      @ct.timeit_quiet
      def build_Pk():
         ik = np.fft.fftfreq(Ns, d=L/Ns)*2*np.pi
         k = ne.evaluate('sqrt(kx**2 + ky**2 + kz**2)', {'kx':ik[:,None,None], 'ky':ik[None,:,None], 'kz':ik[None,None,:(Ns/2+1)]})
         return Pk.compute(k)*L**3

      print np.where(np.isnan(density_out))[0].size
      Nz = np.count_nonzero(cond)
      amplitude = np.sqrt(build_Pk()[cond]/2) if Pk is not None else (1.0/np.sqrt(2))
      density_out.real[cond] = np.random.randn(Nz) * amplitude
      density_out.imag[cond] = np.random.randn(Nz) * amplitude
      print np.where(np.isnan(density_out))[0].size
    
      # Now we have to fix the Nyquist plane
      hNs = Ns/2
      nyquist = density_out[:, :, hNs]
      Nplane = nyquist.size
      nyquist.flat[:Nplane/2] = np.sqrt(2.0)*nyquist.flat[Nplane:Nplane/2:-1].conj()


  return density_out

@ct.timeit_quiet
def run_generation(input_borg, a_borg, a_ic, cosmo, supersample=1, supergenerate=1, do_lpt2=True, shiftPixel=False, psi_instead=False, needvel=True, func='HU_WIGGLES'):
    """ Generate particles and velocities from a BORG snapshot. Returns a tuple of 
    (positions,velocities,N,BoxSize,scale_factor)."""

    borg_vol = ct.read_borg_vol(input_borg)
    N = borg_vol.density.shape[0]

    cgrowth = CosmoGrowth(**cosmo)

    density, L = ba.half_pixel_shift(borg_vol, doshift=shiftPixel)


    # Compute LPT scaling coefficient
    D1 = cgrowth.D(a_ic)
    D1_0 = D1/cgrowth.D(a_borg)
    Dborg = cgrowth.D(a_borg)/cgrowth.D(1.0)
    print "D1_0=%lg" % D1_0

    if supergenerate>1:
      print("Doing supergeneration (factor=%d)" % supergenerate)
      p = ct.CosmologyPower(**cosmo)
      p.setFunction(func)
      p.normalize(cosmo['SIGMA8']*Dborg)
      density = do_supergenerate(density,mulfac=supergenerate,Pk=p,L=L,h=cosmo['h'])

    lpt = LagrangianPerturbation(-density, L, fourier=True, supersample=supersample)

    # Generate grid
    posq = gen_posgrid(N*supersample, L)
    vel= []
    posx = []
    
    velmul = cgrowth.compute_velmul(a_ic) if not psi_instead else 1
    
    D2 = -3./7 * D1_0**2

    if do_lpt2:
      psi2 = lpt.lpt2('all')
    for j in xrange(3):
        # Generate psi_j (displacement along j)
        print("LPT1 axis=%d" % j)
        psi = D1_0*lpt.lpt1(j)
        psi = psi.reshape((psi.size,))
        if do_lpt2:
          print("LPT2 axis=%d" % j)
          psi += D2 * psi2[j].reshape((psi2[j].size,))
        # Generate posx
        posx.append(((posq[j] + psi)%L).astype(np.float32))
        # Generate vel
        if needvel:
          vel.append((psi*velmul).astype(np.float32)) 

    print("velmul=%lg" % (cosmo['h']*velmul))

    lpt.cube.dhat = lpt.dhat
    density = lpt.cube.irfft()
    density *= (cgrowth.D(1)/cgrowth.D(a_borg))

    return posx,vel,density,N*supergenerate*supersample,L,a_ic,cosmo


@ct.timeit_quiet
def whitify(density, L, cosmo, supergenerate=1, zero_fill=False, func='HU_WIGGLES'):

  N = density.shape[0]
  p = ct.CosmologyPower(**cosmo)
  p.setFunction(func)
  p.normalize(cosmo['SIGMA8'])

  @ct.timeit_quiet
  def build_Pk():
    ik = np.fft.fftfreq(N, d=L/N)*2*np.pi
    k = np.sqrt(ik[:,None,None]**2 + ik[None,:,None]**2 + ik[None,None,:(N/2+1)]**2)
    return p.compute(k)*L**3
    
  Pk = build_Pk()
  Pk[0,0,0]=1

  cube = CubeFT(L, N)
  cube.density = density
  density_hat = cube.rfft()
  density_hat /= np.sqrt(Pk)

  Ns = N*supergenerate

  density_hat_super = do_supergenerate(density_hat, mulfac=supergenerate)

  cube = CubeFT(L, Ns)
  cube.dhat = density_hat_super
  return np.fft.irfftn(density_hat_super)*Ns**1.5


    
def write_icfiles(*generated_ic, **kwargs):
  """Write the initial conditions from the tuple returned by run_generation"""
  
  supergenerate=kwargs.get('supergenerate', 1)
  zero_fill=kwargs.get('zero_fill', False)
  posx,vel,density,N,L,a_ic,cosmo = generated_ic

  ct.simpleWriteGadget("Data/borg.gad", posx, velocities=vel, boxsize=L, Hubble=cosmo['h'], Omega_M=cosmo['omega_M_0'], time=a_ic) 
  for i,c in enumerate(["z","y","x"]):
    ct.writeGrafic("Data/ic_velc%s" % c, vel[i].reshape((N,N,N)), L, a_ic, **cosmo)

  ct.writeGrafic("Data/ic_deltab", density, L, a_ic, **cosmo)
  
  ct.writeWhitePhase("Data/white.dat", whitify(density, L, cosmo, supergenerate=supergenerate,zero_fill=zero_fill))
  
  with file("Data/white_params", mode="w") as f:
    f.write("4\n%lg, %lg, %lg\n" % (cosmo['omega_M_0'], cosmo['omega_lambda_0'], 100*cosmo['h']))
    f.write("%lg\n%lg\n-%lg\n0,0\n" % (cosmo['omega_B_0'],cosmo['ns'],cosmo['SIGMA8']))
    f.write("-%lg\n1\n0\n\n\n\n\n" % L)
    f.write("2\n\n0\nwhite.dat\n0\npadding_white.dat\n")

