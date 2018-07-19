import numpy as np

Mpc=3.08e22
rhoc = 1.8864883524081933e-26  # m^(-3)
sigmaT = 6.6524e-29
mp = 1.6726e-27
lightspeed = 299792458.
v_unit = 1e3  # Unit of 1 km/s
T_cmb=2.725
h = 0.71
Y = 0.245  #The Helium abundance
Omega_matter = 0.26
Omega_baryon=0.0445

G=6.67e-11
MassSun=2e30
frac_electron = 1.0 # Hmmmm
frac_gas_galaxy = 0.14
mu = 1/(1-0.5*Y)

tmpp_cat={'Msun':3.29,
    'alpha':-0.7286211634758224,
    'Mstar':-23.172904033796893,
    'PhiStar':0.0113246633636846,
    'lbar':393109973.22508669}

baryon_fraction = Omega_baryon / Omega_matter

ksz_normalization = -T_cmb*sigmaT*v_unit/(lightspeed*mu*mp) * baryon_fraction
assert ksz_normalization < 0
rho_mean_matter = Omega_matter * (3*(100e3/Mpc)**2/(8*np.pi*G)) 
Lbar = tmpp_cat['lbar'] / Mpc**3
M_over_L_galaxy = rho_mean_matter / Lbar


del np
