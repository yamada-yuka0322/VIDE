import numpy as np
import numexpr as ne
from .constants import *

# -----------------------------------------------------------------------------
# Generic profile generator
# -----------------------------------------------------------------------------

class KSZ_Profile(object):
  R_star= 0.0 # 15 kpc/h 
  L_gal0 = 10**(0.4*(tmpp_cat['Msun']-tmpp_cat['Mstar']))

  def __init__(self,sculpt):
    """Base class for KSZ profiles

    Arguments:
      sculpt (float): If negative, do not sculpt. If positive, there will be a 2d 
                      suppression of the profile with a radius given by sculpt (in arcmins).
    """
    self.sculpt = sculpt * np.pi/180/60.
    self.rGalaxy = 1.0

  def evaluate_profile(self, r):
    raise NotImplementedError("Abstract function")

  def projected_profile(self, cos_theta,angularDistance):

    idx_base = idx = np.where(cos_theta > 0)[0]
    tan_theta_2 = 1/(cos_theta[idx]**2) - 1
    tan_theta_2_max = (self.rGalaxy/angularDistance)**2
    tan_theta_2_min = (self.R_star/angularDistance)**2

    idx0 = np.where((tan_theta_2 < tan_theta_2_max))
    idx = idx_base[idx0]
    tan_theta_2 = tan_theta_2[idx0]
    tan_theta = np.sqrt(tan_theta_2)

    r = (tan_theta*angularDistance)

    m,idx_mask = self.evaluate_profile(r)
    idx_mask = idx[idx_mask]

    idx_mask = np.append(idx_mask,idx[np.where(tan_theta_2<tan_theta_2_min)[0]])
    if tan_theta_2.size > 0:
      idx_mask = np.append(idx_mask,idx[tan_theta_2.argmin()])

    if self.sculpt > 0:
      theta = np.arctan(tan_theta)
      cond = theta < self.sculpt
      m[cond] *= (theta[cond]/self.sculpt)**2

    return idx,idx_mask,m
     

# -----------------------------------------------------------------------------
# Isothermal profile generator
# -----------------------------------------------------------------------------

class KSZ_Isothermal(KSZ_Profile):
  sigma_FP=160e3  #m/s
  R_innergal = 0.030

  def __init__(self, Lgal, x, y=0.0, sculpt=-1):
    """Support for Isothermal profile

    Arguments:
      Lgal (float): Galaxy luminosity in solar units
      x (float): extent of halo in virial radius units

    Keyword arguments:
      y (float): Inner part where there is no halo
      sculpt (float): If negative, do not sculpt. If positive, there will be a 2d 
                      suppression of the profile with a radius given by sculpt (in arcmins).
    """
    
    super(KSZ_Isothermal,self).__init__(sculpt)
    
    self.R_gal = 0.226 * x
    self.R_innergal *= y

    self.rho0 = self.sigma_FP**2/(2*np.pi*G) # * (Lgal/L_gal0)**(2./3)
    self.rGalaxy = self.R_gal*(Lgal/self.L_gal0)**(1./3)
    self.rInnerGalaxy = self.R_innergal*(Lgal/self.L_gal0)**(1./3)
    self._prepare()
    
  def _prepare(self):
    pass

  def evaluate_profile(self,r):
    rho0, rGalaxy, rInner = self.rho0, self.rGalaxy, self.rInnerGalaxy
    
    D = {'rho0':rho0, 'rGalaxy':rGalaxy, 'rInner': rInner, 'Mpc':Mpc }
    
    Q = np.zeros(r.size)

    cond = (r<=1e-10)
    Q[cond] = rho0*2/Mpc * (rGalaxy-rInner)/(rGalaxy*rInner)

    cond = (r>0)*(r <= rInner)
    D['r'] = r[cond]
    Q[cond] = ne.evaluate('rho0*2/(Mpc*r) * (arctan(sqrt( (rGalaxy/r)**2 -1 ))  - arctan(sqrt( (rInner/r)**2 - 1 )))',
                local_dict=D)
 
    cond = (r > rInner)*(r <= rGalaxy)
    D['r'] = r[cond]
    Q[cond] = ne.evaluate('rho0*2/(Mpc*r) * arctan(sqrt( (rGalaxy/r)**2 -1 ))',
                local_dict=D)

    return Q,[] #np.where(r<rInner)[0]


# -----------------------------------------------------------------------------
# NFW profile generator
# -----------------------------------------------------------------------------

class KSZ_NFW(KSZ_Profile):
  """ Support for NFW profile
"""

  def __init__(self,x,y=0.0):
    from numpy import log, pi

    if 'pre_nfw' not in self:
      self._prepare()

    kiso = KSZ_Isothermal(x,y)
    r_is = kiso.rGalaxy
    rho_is = kiso.rho0
    r_inner = kiso.rInnerGalaxy

    self.Mgal = rho_is*4*pi*(r_is/args.x)*Mpc #Lgal*M_over_L_galaxy
    self.Rvir = r_is/x

    cs = self._get_concentration(Mgal)

    self.rs = Rvir/cs
    b = (log(1.+cs)-cs/(1.+cs))
    
    self.rho_s = Mgal/(4*pi*b*(rs*Mpc)**3)


  def _prepare(self, _x_min=1e-4, _x_max=1e4):
    from scipy.integrate import quad
    from numpy import sqrt, log10
    from scipy.interpolate import interp1d

    lmin = log10(x_min)
    lmax = log10(x_max)

    x = 10**(np.arange(100)*(lmax-lmin)/100.+lmin)
    profile = np.empty(x.size)
    
    nu_tilde = lambda u: (1/(u*(1+u)**2))

    for i in range(x.size):
        if x[i] < args.x:
            profile[i] = 2*quad(lambda y: (nu_tilde(sqrt(x[i]**2+y**2))), 0, np.sqrt((args.x)**2-x[i]**2))[0]
        else:
            profile[i] = 0

    # Insert the interpolator into the class definition
    KSZ_NFW.pre_nfw = self.pre_nfw = interp1d(x,prof)

  def _get_concentration(self, Mvir):
  	from numpy import exp, log

  	return exp(0.971 - 0.094*log(Mvir/(1e12*MassSun)))

  def evaluate_profile(self,r):
    cs = self._get_concentration(self.Mvir)
    rs = self.Rvir/cs
    
    return self.rho_s*rs*Mpc*self.pre_nfw(r/rs),np.array([],dtype=int)


