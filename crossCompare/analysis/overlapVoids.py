#!/usr/bin/env python
#+
#   VIDE -- Void IDentification and Examination -- ./crossCompare/analysis/mergerTree.py
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

# computes the overlap between two void catalogs

from void_python_tools.backend import *
from void_python_tools.plotting import *
import imp
import pickle
import os
import matplotlib.pyplot as plt
import numpy as np
import argparse

# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='Analyze.')
parser.add_argument('--parm', dest='parm', default='datasetsToAnalyze.py',
                    help='path to parameter file')
args = parser.parse_args()

# ------------------------------------------------------------------------------

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

outFileName = outputDir + "/" + "voidOverlap" #+ ".dat"

with open(workDir+baseSampleDir+"/sample_info.dat", 'rb') as input:
  baseSample = pickle.load(input)

for (iSample, sampleDir) in enumerate(sampleDirList):

  with open(workDir+sampleDir+"/sample_info.dat", 'rb') as input:
    sample = pickle.load(input)

  print " Working with", sample.fullName, "...",
  sys.stdout.flush()

  sampleName = sample.fullName

  binPath = CTOOLS_PATH+"/analysis/voidOverlap"
  logFile = logDir+"/mergertree_"+baseSample.fullName+"_"+sampleName+".out"
  stepOutputFileName = outFileName + "_" + baseSample.fullName + "_" + \
                       sampleName + "_"

  launchVoidOverlap(baseSample, sample, workDir+baseSampleDir, 
                    workDir+sampleDir, binPath, 
                    thisDataPortion="central", logFile=logFile,
                    continueRun=False, workDir=workDir,
                    outputFile=stepOutputFileName,
                    #matchMethod="useID")
                    matchMethod="prox")

print " Done!"
