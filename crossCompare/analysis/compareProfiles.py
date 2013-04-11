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

if not os.access(logDir, os.F_OK):
  os.makedirs(logDir)

mergerFileBase = dataDir + "/" + mergerNameBase

# get list of base voids
with open(workDir+baseSampleDir+"/sample_info.dat", 'rb') as input:
  baseSample = pickle.load(input)
  baseSampleName = baseSample.fullName
  baseVoidList = np.loadtxt(workDir+baseSampleDir+"/untrimmed_centers_central_"+\
                    baseSampleName+".out")

for stack in baseSample.stacks:
  print " Stack:", stack.rMin    

  accepted = (baseVoidList[:,4] > stack.rMin) & (baseVoidList[:,4] < stack.rMax)
  baseIDList = baseVoidList[:,7][accepted]

  for (iSample, sampleDir) in enumerate(sampleDirList):
    with open(workDir+sampleDir+"/sample_info.dat", 'rb') as input:
      sample = pickle.load(input)

    print "  Working with", sample.fullName, "..."
    sys.stdout.flush()
    sampleName = sample.fullName

    # get list of appropriate voids
    if sample.fullName == baseSample.fullName:
      idList = baseIDList
    else:
      matchList = np.loadtxt(mergerFileBase+"_"+baseSampleName+"_"+sampleName+\
                            "_summary.out")
      idList = []
      for i,testID in enumerate(matchList[:,8]):
        if np.any(testID == baseIDList):
          idList.append(matchList[i,0])
      idList = np.array(idList)
      
    idList = idList.astype(int)

    print "    Found", len(idList), "voids to work with"
 
    voidBaseDir = workDir+"/"+sampleDir+"stacks"
        
    runSuffix = getStackSuffix(stack.zMin, stack.zMax, stack.rMin, 
                               stack.rMax, dataPortion, 
                               customLine="selected")
    voidDir = voidBaseDir+"_"+runSuffix

    if not os.access(voidDir,os.F_OK): os.makedirs(voidDir)

    if len(idList) == 0:
      print "   No voids here anyway, skipping..."
      continue
  
    print "    Stacking voids...",
    sys.stdout.flush()
        
    STACK_PATH = CTOOLS_PATH+"/stacking/stackVoidsZero"

    launchStack(sample, stack, STACK_PATH, 
                   thisDataPortion=dataPortion,
                   logDir=logDir, voidDir=voidDir, 
                   zobovDir=workDir+"/"+sampleDir,
                   freshStack=True, INCOHERENT=False,
                   ranSeed=101010, summaryFile=None, 
                   continueRun=False,
                   dataType=sample.dataType,
                   idList=idList, rescaleOverride="")

    print "    Profiling stacks...",
    sys.stdout.flush()
    logFile = logDir+"/profile_"+sampleName+"_"+runSuffix+".out"
    launchProfile(sample, stack, voidDir=voidDir, 
                  logFile=logFile, continueRun=False)

print " Done!"
