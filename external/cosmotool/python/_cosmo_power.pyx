from libcpp cimport bool
from libcpp cimport string as cppstring
import numpy as np
cimport numpy as np
from cpython cimport PyObject, Py_INCREF
cimport cython

np.import_array()

cdef extern from "cosmopower.hpp" namespace "CosmoTool":

  cdef enum CosmoFunction "CosmoTool::CosmoPower::CosmoFunction":
      POWER_EFSTATHIOU "CosmoTool::CosmoPower::POWER_EFSTATHIOU",
      HU_WIGGLES "CosmoTool::CosmoPower::HU_WIGGLES",
      HU_BARYON "CosmoTool::CosmoPower::HU_BARYON",
      OLD_POWERSPECTRUM,
      POWER_BARDEEN "CosmoTool::CosmoPower::POWER_BARDEEN",
      POWER_SUGIYAMA "CosmoTool::CosmoPower::POWER_SUGIYAMA",
      POWER_BDM,
      POWER_TEST,
      HU_WIGGLES_ORIGINAL "CosmoTool::CosmoPower::HU_WIGGLES_ORIGINAL"

  cdef cppclass CosmoPower:
    double n
    double K0
    double V_LG_CMB
    
    double CMB_VECTOR[3]
    double h
    double SIGMA8
    double OMEGA_B
    double OMEGA_C
    double omega_B
    double omega_C
    double Theta_27
    double OMEGA_0
    double Omega
    double beta
    double OmegaEff
    double Gamma0
    double normPower
    
        
    CosmoPower()
    void setFunction(CosmoFunction)
    void updateCosmology()
    void updatePhysicalCosmology()
    void normalize(double,double)
    void setNormalization(double)
    double power(double)
 
cdef class CosmologyPower:
  """CosmologyPower(**cosmo)
  
  CosmologyPower manages and compute power spectra computation according to different
  approximation given in the litterature.
  
  Keyword arguments:
    omega_B_0 (float): relative baryon density
    omega_M_0 (float): relative matter density
    h (float): Hubble constant relative to 100 km/s/Mpc
    ns (float): power law of the large scale inflation spectrum
  """

  cdef CosmoPower power

  def __init__(self,**cosmo):    
    self.power = CosmoPower()
    self.power.OMEGA_B = cosmo['omega_B_0']
    self.power.OMEGA_C = cosmo['omega_M_0']-cosmo['omega_B_0']
    self.power.h = cosmo['h']
    if 'ns' in cosmo:
      self.power.n = cosmo['ns'] 
    
    assert self.power.OMEGA_C > 0
    
    self.power.updateCosmology()

  def setNormalization(self,A):
    self.power.setNormalization(A)
    
  def normalize(self,s8,k_min=-1,k_max=-1):
    """normalize(self, sigma8)
    
    Compute the normalization of the power spectrum using sigma8.
    
    Arguments:
      sigma8 (float): standard deviation of density field smoothed at 8 Mpc/h
    """
    self.power.SIGMA8 = s8
    self.power.normalize(k_min, k_max)
    
    
  def setFunction(self,funcname):
    """setFunction(self, funcname)
    
    Choose an approximation to use for the computation of the power spectrum
    
    Arguments:
      funcname (str): the name of the approximation. It can be either
                      EFSTATHIOU, HU_WIGGLES, HU_BARYON, BARDEEN or SUGIYAMA.
    """
    cdef CosmoFunction f
    
    f = POWER_EFSTATHIOU
  
    if funcname=='EFSTATHIOU':
      f = POWER_EFSTATHIOU
    elif funcname=='HU_WIGGLES':
      f = HU_WIGGLES
    elif funcname=='HU_BARYON':
      f = HU_BARYON
    elif funcname=='BARDEEN':
      f = POWER_BARDEEN
    elif funcname=='SUGIYAMA':
      f = POWER_SUGIYAMA
    elif funcname=='HU_WIGGLES_ORIGINAL':
      f = HU_WIGGLES_ORIGINAL
    else:
      raise ValueError("Unknown function name " + funcname)

    self.power.setFunction(f)
    
  cdef double _compute(self, double k):
    k *= self.power.h
    return self.power.power(k) * self.power.h**3
    
  def compute(self, k):
    """compute(self, k)
    
    Compute the power spectrum for mode which length k.
    
    Arguments:
      k (float): Mode for which to evaluate the power spectrum.
                 It can be a scalar or a numpy array. 
                 The units must be in 'h Mpc^{-1}'.
    
    Returns:
      a scalar or a numpy array depending on the type of the k argument 
    """
    
    cdef np.ndarray out
    cdef double kval
    cdef tuple i
    
    if isinstance(k, np.ndarray):
      out = np.empty(k.shape, dtype=np.float64)
      for i,kval in np.ndenumerate(k):
        out[i] = self._compute(kval)
      return out
    else:
      return self._compute(k)
      
