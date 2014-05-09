#+
#   VIDE -- Void IDentification and Examination -- ./python_tools/void_python_tools/plotting/plotTools.py
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
__all__=['buildProfile','fitHSWProfile','getHSWProfile',]

from void_python_tools.backend.classes import *
from plotDefs import *
import numpy as np
import os
import void_python_tools.apTools as vp
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

def HamausProfile(r, rs, dc):
  alpha = -2.0*rs + 4.0
  if rs < 0.91:
    beta = 17.5*rs - 6.5
  else:
    beta = -9.8*rs + 18.4
  return dc * (1 - (r/rs)**alpha) / (1+ (r)**beta) + 1

# -----------------------------------------------------------------------------
def buildProfile(catalog, rMin, rMax):

# builds a stacked radial density profile from the given catalog
#   catalog: void catalog
#   rMin: minimum void radius, in Mpc/h
#   rMax: maximum void radius, in Mpc/h
#
# returns:
#   binCenters: array of radii in binned profile
#   stackedProfile: the stacked density profile
#   sigmas: 1-sigma uncertainty in each bin

  rMaxProfile = rMin*3 + 2
  periodicLine = getPeriodic(catalog.sampleInfo)

  print "  Building particle tree..."
  partTree = getPartTree(catalog)

  print "  Selecting voids to stack..."
  accepted = (catalog.voids[:].radius > rMin) & (catalog.voids[:].radius < rMax)
  voidsToStack = catalog.voids[accepted]

  print "  Stacking voids..."
  allProfiles = []
  for void in voidsToStack:
    center = void.barycenter
    
    localPart = catalog.partData[ getBall(partTree, center, rMaxProfile) ]
    shiftedPart = shiftPart(localPart, center, periodicLine, catalog.ranges)

    dist = np.sqrt(np.sum(shiftedPart[:,:]**2, axis=1))
    thisProfile, radii = np.histogram(dist, bins=nBins, range=(0,rMaxProfile))
    deltaV = 4*np.pi/3*(radii[1:]**3-radii[0:(radii.size-1)]**3)
    thisProfile = np.float32(thisProfile)
    thisProfile /= deltaV
    thisProfile /= catalog.volNorm
    allProfiles.append(thisProfile)
    binCenters = 0.5*(radii[1:] + radii[:-1])

  stackedProfile = np.std(allProfiles, axis=0) / np.sqrt(nVoids)
  sigmas = np.std(allProfiles, axis=0) / np.sqrt(nVoids)
  
  return binCenters, stackedProfile, sigmas


# -----------------------------------------------------------------------------
def fitHSWProfile(radii, densities, sigmas):

# fits the given density profile to the HSW function
#   radii: array of radii in r/rV
#   densities: array of densities
#   sigmas: array of uncertainties
# 
# returns:
#   popt: best-fit values of dc and rs
#   pcov: covariance matrix

  popt, pcov = curve_fit(HamausProfile, radii, densities,
                         sigma=sigmas)
                         maxfev=10000, xtol=5.e-3,
                         p0=[1.0,-1.0])

  return popt, pcov


# -----------------------------------------------------------------------------
def getHSWProfile(den, radius):

# returns the HSW profile for the given sample density and void size
# (interpolated from best-fit values)
#  density: density of sample
#  radius: void size in Mpc/h

# returns:
#   binCenters: array of radii in binned profile
#   stackedProfile: the density profile

  sample = catalog.sampleInfo
  data = catalog.voids[:].radius

