#+
#   VIDE -- Void IDentification and Examination -- ./python_tools/void_python_tools/apTools/chi2/cosmologyTools.py
#   Copyright (C) 2010-2013 Guilhem Lavaux
#   Copyright (C) 2011-2013 P. M. Sutter
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; version 2 of the License.
# 
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program; if not, write to the Free Software Foundation, Inc.,
#   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#+
# a suite of functions to compute expansion rates, angular diameter 
# distances, and expected void stretching

import numpy as np
import scipy.integrate as integrate
from void_python_tools.backend import *

__all__=['expansion', 'angularDiameter', 'expectedStretch', 'aveStretch', 'aveExpansion', 'aveStretchCone', 'aveWeightedStretch']

# returns 1/E(z) for the given cosmology
def expansion(z, Om = 0.27, Ot = 1.0, w0 = -1.0, wa = 0.0):
  ez = Om * (1+z)**3 + (Ot-Om)# * (1+z)**(3.+3*wz)
  #ez = Om * (1+z)**3 + (Ot-Om)# * integrade.quad(eosDE, 0.0, z, args=(w0,wa))[0]
  ez = 1./np.sqrt(ez)
  return ez

# returns DE value at redshift z
def eosDE(z, w0 = -1.0, wa = 0.0):
  return w0 + wa*z/(1+z)
  
# returns D_A(z) for the given cosmology
def angularDiameter(z, Om = 0.27, Ot = 1.0, w0 = -1.0, wa = 0.0):
  da = integrate.quad(expansion, 0.0, z, args=(Om, Ot, w0, wa))[0]
  return da

  
# returns expected void stretch for the given cosmology
def expectedStretch(z, Om = 0.27, Ot = 1.0, w0 = -1.0, wa = 0.0):
  ez = 1./expansion(z, Om=Om, Ot=Ot, w0=w0, wa=wa)
  da = angularDiameter(z, Om=Om, Ot=Ot, w0=w0, wa=wa)
  return ez*da/z

# -----------------------------------------------------------------------------
# returns average expected void stretch for a given redshift range 
#   assuming a cone 
def aveStretchCone(zStart, zEnd, skyFrac = 0.19, Om = 0.27, Ot = 1.0, 
                   w0 = -1.0, wa = 0.0):
  #print "assuming observation!", skyFrac
  if zStart == 0.0: zStart = 1.e-6

  h1 = zStart 
  h2 = zEnd

  r1 = skyFrac * 4* np.pi * zStart**2
  r2 = skyFrac * 4 * np.pi * zEnd**2

  # surface area of a slice within a cone
  def coneSlice(x, h, r):
    return np.pi * (r/h*x)**2

  def coneFunc(z, h, r, Om = 0.27, Ot = 1.0, w0 = -1.0, wa = 0.0):
    return coneSlice(z, h, r) * expectedStretch(z, Om, Ot, w0, wa)

  aveHigh = integrate.quad(coneFunc, 0.0, zEnd, args=(h2, r2, Om, Ot, w0, wa), full_output=1)[0]
  aveLow = integrate.quad(coneFunc, 0.0, zStart, args=(h1, r1, Om, Ot, w0, wa), full_output=1)[0]
  volumeHigh = integrate.quad(coneSlice, 0.0, zEnd, args=(h2, r2))[0]
  volumeLow = integrate.quad(coneSlice, 0.0, zStart, args=(h1, r1))[0]

  return (aveHigh-aveLow)/(volumeHigh-volumeLow)

# returns average expected void stretch for a given redshift range
def aveStretch(sample, zStart, zEnd, 
               Om = 0.27, Ot = 1.0, w0 = -1.0, wa = 0.0):
  if zStart == 0.0: zStart = 1.e-6

  if sample.dataType == "observation":
    stretch =  aveStretchCone(zStart, zEnd, 
                          skyFrac=sample.skyFraction, Om=Om, Ot=Ot, 
                          w0=w0, wa=wa)
  else:
    ave = integrate.quad(expectedStretch, zStart, zEnd, 
                         args=(Om, Ot, w0, wa))[0]
    ave /= (zEnd-zStart)
    stretch = ave

  # if in comoving space, calculate stretch for fiducial cosmology
  # and take relative amount
  if not sample.useLightCone or sample.useComoving:
    if sample.dataType == "observation":
      stretchFid =  aveStretchCone(zStart, zEnd, 
                          skyFrac=sample.skyFraction, Om=sample.omegaM, Ot=Ot, 
                          w0=w0, wa=wa)
    else:
      ave = integrate.quad(expectedStretch, zStart, zEnd, 
                           args=(sample.omegaM, Ot, w0, wa))[0]
      ave /= (zEnd-zStart)
      stretchFid = ave

  stretch = stretchFid/stretch

  return stretch


# -----------------------------------------------------------------------------
# returns average expected void stretch for a given redshift range 
#   weighted by an actual void distribution
def aveWeightedStretch(zStart, zEnd, skyFrac = 0.19, Om = 0.27, Ot = 1.0, 
                   w0 = -1.0, wa = 0.0, dist=None, bins=None):
  if zStart == 0.0: zStart = 1.e-6

  def weightedSlice(x):
    return np.interp(x, bins[:-1], dist)
 
  def weightedFunc(z, Om = 0.27, Ot = 1.0, w0 = -1.0, wa = 0.0):
    return expectedStretch(z, Om, Ot, w0, wa) *\
           weightedSlice(z) 

  ave = integrate.quad(weightedFunc, zStart, zEnd, args=(Om, Ot, w0, wa),
                       full_output=1)[0]

  volume = integrate.quad(weightedSlice, zStart, zEnd, full_output=1)[0]
   
  return ave/volume


# -----------------------------------------------------------------------------
# returns average expected expansion for a given redshift range
def aveExpansion(zStart, zEnd, Om = 0.27, Ot = 1.0, w0 = -1.0, wa = 0.0):
  if zStart == 0.0: zStart = 1.e-6
  ave = integrate.quad(expansion, zStart, zEnd, args=(Om, Ot, w0, wa))[0]
  ave = (zEnd-zStart)/ave
  return ave


