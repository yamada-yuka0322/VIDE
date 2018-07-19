import numexpr as ne
import multiprocessing
import pyfftw
import weakref
import numpy as np
import cosmolopy as cpy
import cosmotool as ct

class CubeFT(object):
  def __init__(self, L, N, max_cpu=-1):
  
      self.N = N
      self.align = pyfftw.simd_alignment
      self.L = L
      self.max_cpu = multiprocessing.cpu_count() if max_cpu < 0 else max_cpu
      self._dhat = pyfftw.n_byte_align_empty((self.N,self.N,self.N/2+1), self.align, dtype='complex64')
      self._density = pyfftw.n_byte_align_empty((self.N,self.N,self.N), self.align, dtype='float32')
      self._irfft = pyfftw.FFTW(self._dhat, self._density, axes=(0,1,2), direction='FFTW_BACKWARD', threads=self.max_cpu, normalize_idft=False)
      self._rfft = pyfftw.FFTW(self._density, self._dhat, axes=(0,1,2), threads=self.max_cpu, normalize_idft=False)

  def rfft(self):
      return ne.evaluate('c*a', local_dict={'c':self._rfft(normalise_idft=False),'a':(self.L/self.N)**3})
      
  def irfft(self):
      return ne.evaluate('c*a', local_dict={'c':self._irfft(normalise_idft=False),'a':(1/self.L)**3})
   
  def get_dhat(self):
      return self._dhat
  def set_dhat(self, in_dhat):
      self._dhat[:] = in_dhat 
  dhat = property(get_dhat, set_dhat, None)
      
  def get_density(self):
      return self._density      
  def set_density(self, d):
      self._density[:] = d 
  density = property(get_density, set_density, None)
      

class CosmoGrowth(object):

  def __init__(self, **cosmo):
    self.cosmo = cosmo

  def D(self, a):
    return cpy.perturbation.fgrowth(1/a-1, self.cosmo['omega_M_0'], unnormed=True)

  def compute_E(self, a):
    om = self.cosmo['omega_M_0']
    ol = self.cosmo['omega_lambda_0']
    ok = self.cosmo['omega_k_0']

    E = np.sqrt(om/a**3 + ol + ok/a**2)
  
    H2 = -3*om/a**4  - 2*ok/a**3

    Eprime = 0.5*H2/E

    return E,Eprime

  def Ddot(self, a):
    E,Eprime = self.compute_E(a)
    D = self.D(a)
    Ddot_D = Eprime/E + 2.5 * self.cosmo['omega_M_0']/(a**3*E**2*D)
    Ddot_D *= a
    return Ddot_D

  def compute_velmul(self, a):
    E,_ = self.compute_E(a)
    velmul = self.Ddot(a)
    velmul *= 100 * a * E
    return velmul





class LagrangianPerturbation(object):

    def __init__(self,density,L, fourier=False, supersample=1, max_cpu=-1):
    
        self.L = L
        self.N = density.shape[0]
        
        self.max_cpu = max_cpu
        self.cube = CubeFT(self.L, self.N, max_cpu=max_cpu)
        
        if not fourier:
          self.cube.density = density
          self.dhat = self.cube.rfft().copy()
        else:
          self.dhat = density.copy()

        if supersample > 1:
          self.upgrade_sampling(supersample)
        self.ik = np.fft.fftfreq(self.N, d=L/self.N)*2*np.pi
        self._kx = self.ik[:,None,None]
        self._ky = self.ik[None,:,None]
        self._kz = self.ik[None,None,:(self.N/2+1)]
        self.cache = {}#weakref.WeakValueDictionary()

    @ct.timeit_quiet
    def upgrade_sampling(self, supersample):
        N2 = self.N * supersample
        N = self.N
        dhat_new = np.zeros((N2, N2, N2/2+1), dtype=np.complex128)

        hN = N/2
        dhat_new[:hN, :hN, :hN+1] = self.dhat[:hN, :hN, :]
        dhat_new[:hN, (N2-hN):N2, :hN+1] = self.dhat[:hN, hN:, :]
        dhat_new[(N2-hN):N2, (N2-hN):N2, :hN+1] = self.dhat[hN:, hN:, :]
        dhat_new[(N2-hN):N2, :hN, :hN+1] = self.dhat[hN:, :hN, :]

        self.dhat = dhat_new
        self.N = N2
        self.cube = CubeFT(self.L, self.N, max_cpu=self.max_cpu)
        
    @ct.timeit_quiet
    def _gradient(self, phi, direction):
        if direction == 'all':
          dirs = [0,1,2]
          copy = True
        else:
          dirs = [direction]
          copy = False
        ret=[]
        for dir in dirs:
          ne.evaluate('phi_hat * i * kv / (kx**2 + ky**2 + kz**2)', out=self.cube.dhat,
              local_dict={'i':-1j, 'phi_hat':phi, 'kv':self._kdir(dir),
                          'kx':self._kx, 'ky':self._ky, 'kz':self._kz},casting='unsafe')
#        self.cube.dhat = self._kdir(direction)*1j*phi 
          self.cube.dhat[0,0,0] = 0
          x = self.cube.irfft()
          ret.append(x.copy() if copy else x)
        return ret[0] if len(ret)==1 else ret
            
    @ct.timeit_quiet
    def lpt1(self, direction=0):
        return self._gradient(self.dhat, direction)
        
    def new_shape(self,direction, q=3, half=False):
        N0 = (self.N/2+1) if half else self.N
        return ((1,)*direction) + (N0,) + ((1,)*(q-1-direction))

    def _kdir(self, direction, q=3):
        if direction != q-1:
            return self.ik.reshape(self.new_shape(direction, q=q))
        else:
            return self.ik[:self.N/2+1].reshape(self.new_shape(direction, q=q, half=True))

    def _get_k2(self, q=3):
        if 'k2' in self.cache:
            return self.cache['k2']
            
        k2 = self._kdir(0, q=q)**2
        for d in xrange(1,q):
            k2 = k2 + self._kdir(d, q=q)**2
            
        self.cache['k2'] = k2
        return k2
        
    def _do_irfft(self, array, copy=True):
        if copy:
          self.cube.dhat = array
        return self.cube.irfft()
        
    def _do_rfft(self, array, copy=True):
        if copy:
          self.cube.density = array
        return self.cube.rfft()

    @ct.timeit_quiet
    def lpt2(self, direction=0):
#        k2 = self._get_k2()
#        k2[0,0,0] = 1

        inv_k2 = ne.evaluate('1/(kx**2+ky**2+kz**2)', {'kx':self._kdir(0),'ky':self._kdir(1),'kz':self._kdir(2)})
        inv_k2[0,0,0]=0
        potgen0 = lambda i: ne.evaluate('kdir**2*d*ik2',out=self.cube.dhat,local_dict={'kdir':self._kdir(i),'d':self.dhat,'ik2':inv_k2}, casting='unsafe' )
        potgen = lambda i,j: ne.evaluate('kdir0*kdir1*d*ik2',out=self.cube.dhat,local_dict={'kdir0':self._kdir(i),'kdir1':self._kdir(j),'d':self.dhat,'ik2':inv_k2}, casting='unsafe' )

        if 'lpt2_potential' not in self.cache:
            print("Rebuilding potential...")
            div_phi2 = np.zeros((self.N,self.N,self.N), dtype=np.float64)
            for j in xrange(3):
                q = self._do_irfft( potgen0(j) ).copy()
                for i in xrange(j+1, 3):
                    with ct.time_block("LPT2 elemental (%d,%d)" %(i,j)):
                      ne.evaluate('div + q * pot', out=div_phi2, 
                            local_dict={'div':div_phi2, 'q':q,'pot':self._do_irfft( potgen0(i), copy=False ) } 
                            )
                      ne.evaluate('div - pot**2',out=div_phi2,
                            local_dict={'div':div_phi2,'pot':self._do_irfft(potgen(i,j), copy=False) }
                            )

            phi2_hat = self._do_rfft(div_phi2)
            #self.cache['lpt2_potential'] = phi2_hat
            del div_phi2
        else:
            phi2_hat = self.cache['lpt2_potential']
            
        return self._gradient(phi2_hat, direction)
