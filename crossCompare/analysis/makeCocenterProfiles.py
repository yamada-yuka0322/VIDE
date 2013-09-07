#!/usr/bin/env python
#+
#   VIDE -- Void IDEntification pipeline -- ./pipeline/apAnalysis.py
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

# takes voids of given radii, computes 1 profiles, 
#   then computes 1d profiles  for higher-resolution catalogs using
#   same positions

# computes radial density profiles centered on baseSample

import imp
import pickle
import os
import numpy as np
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from void_python_tools.backend import *
from util import *

# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='Analyze.')
parser.add_argument('--parm', dest='parm', default='datasetsToAnalyze.py', help='path to parameter file')
parser.add_argument('--show', dest='showPlot', action='store_const',
                   const=True, default=False,
                   help='display the plot (default: just write eps)')
args = parser.parse_args()

# -----------------------------------------------------------------------------
# plot a slice of the density around the void in baseIDList, 
#   with any voids in the slice shown and any voids in baseIDList flagged
def saveProfiles(baseSample, stack, sampleList, profileList, radii,
                 figDir, showPlot, outputDir):

  thisRadius = str(stack.rMin) + "-" + str(stack.rMax)

  plotName = "1dprofile_cocenter_" + baseSample.fullName+"_"+thisRadius

  np.savez(outputDir+"/1dprofile_cocentered_"+plotName+".dat", 
           profileList, radii)

  return

# -----------------------------------------------------------------------------

filename = args.parm
print " Loading parameters from", filename
if not os.access(filename, os.F_OK):
  print "  Cannot find parameter file %s!" % filename
  exit(-1)
parms = imp.load_source("name", filename)
globals().update(vars(parms))

if not os.access(outputDir, os.F_OK):
  os.makedirs(outputDir)

if not os.access(logDir, os.F_OK):
  os.makedirs(logDir)

if not os.access(figDir, os.F_OK):
  os.makedirs(figDir)

# get list of base voids
with open(workDir+baseSampleDir+"/sample_info.dat", 'rb') as input:
  baseSample = pickle.load(input)
  baseSampleName = baseSample.fullName
  baseVoidList = np.loadtxt(workDir+baseSampleDir+"/centers_central_"+\
                    baseSampleName+".out")

sampleList = []
for sampleDir in sampleDirList:
  if compareSampleTag in sampleDir: continue
  with open(workDir+sampleDir+"/sample_info.dat", 'rb') as input:
    sampleList.append(pickle.load(input))

sampleDirList.insert(0,baseSampleDir)
sampleList.insert(0,baseSample)

# pick our void sample
for stack in baseSample.stacks:
  print " Stack:", stack.rMin, "-", stack.rMax
  accepted = (baseVoidList[:,4] > stack.rMin) & (baseVoidList[:,4] < stack.rMax)
  stackVoidList = baseVoidList[accepted]
  print "  We have", len(stackVoidList), "voids here"

  profileList = []
  radii = []

  rMaxProfile = stack.rMin*3 + 2
  if baseSample.profileBinSize == "auto":
    density = 0.5 * 50 / rMaxProfile / 2
  else:
    density = baseSample.profileBinSize
  nBins = rMaxProfile*density

  for (iSample, sampleDir) in enumerate(sampleDirList):
    if compareSampleTag in sampleDir: continue
    sample = sampleList[iSample]
    print "  Working with", sample.fullName, "..."
    sys.stdout.flush()
    sampleName = sample.fullName

    print "   Loading particle data..."
    partData, boxLen, volNorm = loadPart(workDir, sampleDir, sample)

    stackedProfile = np.zeros((nBins))

    print "   Stacking voids..."
    binCenters = [] 
    for void in stackVoidList:
      periodicLine = getPeriodic(sample)
      center = void[0:3]
      shiftedPart = shiftPart(partData, center, periodicLine, boxLen)

      dist = np.sqrt(shiftedPart[:,0]**2 + shiftedPart[:,1]**2 + \
                     shiftedPart[:,2]**2)  
      thisProfile, radii = np.histogram(dist, bins=nBins, range=(0,rMaxProfile))
      deltaV = 4*np.pi/3*(radii[1:]**3-radii[0:(radii.size-1)]**3)   
      thisProfile = np.float32(thisProfile)
      thisProfile /= deltaV
      stackedProfile += thisProfile
      binCenters = 0.5*(radii[1:]+radii[:-1])
      
    stackedProfile /= volNorm
    stackedProfile /= len(stackVoidList)

    profileList.append(stackedProfile)

  # plot these profiles
  print "  Plotting..."
  sys.stdout.flush()
  #binCenters = 0.5*(radii[1:] + radii[:-1])

  saveProfiles(baseSample, stack, sampleList, profileList, binCenters,
               figDir, args.showPlot, outputDir)
    
print " Done!"

