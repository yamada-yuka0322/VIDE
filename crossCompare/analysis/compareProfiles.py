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

# takes a set of voids with given radii, and stacks and plots matches in 
# other catalogs

from void_python_tools.backend import *
import imp
import pickle
import os
import numpy as np
import argparse

# ------------------------------------------------------------------------------

mergerNameBase = "mergerTree"

parser = argparse.ArgumentParser(description='Analyze.')
parser.add_argument('--parm', dest='parm', default='datasetsToAnalyze.py', help='path to parameter file')
args = parser.parse_args()

# ------------------------------------------------------------------------------

filename = args.parm
print " Loading parameters from", filename
if not os.access(filename, os.F_OK):
  print "  Cannot find parameter file %s!" % filename
  exit(-1)
parms = imp.load_source("name", filename)
globals().update(vars(parms))

if not os.access(dataDir, os.F_OK):
  os.makedirs(dataDir)

mergerFileBase = dataDir + "/" + mergerNameBase

baseIDList = []

# get list of base voids
with open(workDir+baseSampleDir+"/sample_info.dat", 'rb') as input:
  baseSample = pickle.load(input)
  baseSampleName = baseSample.fullName
  baseVoidList = np.loadtxt(workDir+baseSampleDir+"/centers_nocut_central_"+\
                    baseSampleName+".out")

for stack in baseSample.stacks:
  print " Stack:", stack.rMin    

  accepted = (baseVoidList[:,4] > stack.rMin) & (baseVoidList[:,4] < stack.rMax)
  baseIDList = baseVoidList[accepted][:,7]

  for (iSample, sampleDir) in enumerate(sampleDirList):
    with open(workDir+sampleDir+"/sample_info.dat", 'rb') as input:
      sample = pickle.load(input)

    print "  Working with", sample.fullName, "...",
    sys.stdout.flush()
    sampleName = sample.fullName

    # get list of appropriate voids
    if sample == baseSample:
      idList = baseIDList
    else:
      idList = 2
      matchList = np.loadtxt(mergerFileBase+"_"+baseSampleName+"_"+sampleName+\
                            "_summary.out")
      accepted = (matchList[:,0] == baseIDList).any()
      idList = matchList[accepted][:,8]

    voidBaseDir = workDir+"/"+sampleDir+"stacks"
        
    runSuffix = getStackSuffix(stack.zMin, stack.zMax, stack.rMin, 
                               stack.rMax, thisDataPortion, 
                               customLine="selected")
    stack.rMin = 0.
    stack.rMax = 1000.

    voidDir = voidBaseDir+"_"+runSuffix

    if not os.access(voidDir,os.F_OK): os.makedirs(voidDir)

    print "    Stacking voids...",
    sys.stdout.flush()
        
    STACK_PATH = CTOOLS_PATH+"/stacking/stackVoidsZero"

    launchStack(sample, stack, STACK_PATH, 
                   thisDataPortion=thisDataPortion,
                   logDir=logDir, voidDir=voidDir, 
                   zobovDir=workDir+"/"+sampleDir,
                   freshStack=freshStack, INCOHERENT=False,
                   ranSeed=101010, summaryFile=None, 
                   continueRun=continueRun,
                   dataType=sample.dataType,
                   idList=idList)
 
    print "   Profiling stacks...",
    sys.stdout.flush()
    logFile = logDir+"/profile_"+sampleName+"_"+runSuffix+".out"
    launchProfile(sample, stack, voidDir=voidDir, 
                  logFile=logFile, continueRun=False)

print " Done!"
