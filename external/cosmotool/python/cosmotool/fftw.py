import pyfftw
import multiprocessing
import numpy as np
import numexpr as ne

class CubeFT(object):
  def __init__(self, L, N, max_cpu=-1, width=32):
  
      if width==32:
        fourier_type='complex64'
        real_type='float32'
      elif width==64:
        fourier_type='complex128'
        real_type='float64'
      else:
        raise ValueError("Invalid bitwidth (must be 32 or 64)")

      self.N = N
      self.align = pyfftw.simd_alignment
      self.L = float(L)
      self.max_cpu = multiprocessing.cpu_count() if max_cpu < 0 else max_cpu
      self._dhat = pyfftw.n_byte_align_empty((self.N,self.N,self.N//2+1), self.align, dtype=fourier_type)
      self._density = pyfftw.n_byte_align_empty((self.N,self.N,self.N), self.align, dtype=real_type)
      self._irfft = pyfftw.FFTW(self._dhat, self._density, axes=(0,1,2), direction='FFTW_BACKWARD', threads=self.max_cpu)#, normalize_idft=False)
      self._rfft = pyfftw.FFTW(self._density, self._dhat, axes=(0,1,2), threads=self.max_cpu) #, normalize_idft=False)

  def rfft(self):
      return ne.evaluate('c*a', out=self._dhat, local_dict={'c':self._rfft(normalise_idft=False),'a':(self.L/self.N)**3}, casting='unsafe')
      
  def irfft(self):
      return ne.evaluate('c*a', out=self._density, local_dict={'c':self._irfft(normalise_idft=False),'a':(1/self.L)**3}, casting='unsafe')
   
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
      

