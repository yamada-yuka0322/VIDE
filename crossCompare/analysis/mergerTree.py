#!/usr/bin/env python

# plots cumulative distributions of number counts

from void_python_tools.backend import *
from void_python_tools.plotting import *
import imp
import pickle
import os
import matplotlib.pyplot as plt
import numpy as np
import argparse

# ------------------------------------------------------------------------------

dataNameBase = "mergerTree"

parser = argparse.ArgumentParser(description='Analyze.')
parser.add_argument('--parmFile', dest='parmFile', default='datasetsToAnalyze.py',
                    help='path to parameter file')
args = parser.parse_args()

# ------------------------------------------------------------------------------

filename = args.parmFile
print " Loading parameters from", filename
if not os.access(filename, os.F_OK):
  print "  Cannot find parameter file %s!" % filename
  exit(-1)
parms = imp.load_source("name", filename)
globals().update(vars(parms))

if not os.access(figDir, os.F_OK):
  os.makedirs(figDir)

outFileName = dataDir + "/" + dataNameBase + ".dat"

for (iSample, sampleDir) in enumerate(sampleDirList):
  if iSample == 0: continue

  with open(workDir+sampleDir+"/sample_info.dat", 'rb') as input:
    sample = pickle.load(input)

  with open(workDir+sampleDirList[0]+"/sample_info.dat", 'rb') as input:
    baseSample = pickle.load(input)

  print " Working with", sample.fullName, "..."

  sampleName = sample.fullName

  binPath = CTOOLS_PATH+"/analysis/voidOverlap"
  logFile = os.getcwd()+"/mergerTree.log"
  stepOutputFileName = os.getcwd()+"/thisStep.out"

  launchVoidOverlap(baseSample, sample, workDir+sampleDirList[0], 
                    workDir+sampleDir, binPath, 
                    thisDataPortion="central", logFile=logFile,
                    continueRun=False, workDir=workDir,
                    outputFile=stepOutputFileName)

  # attach columns to summary file
  if iSample == 1:
    os.system("cp %s %s" % (stepOutputFileName, outFileName))
  else:
    outFile = fp.open("temp.out", "w")
    stepOutFile1 = fp.open(outFileName, "r")
    stepOutFile2 = fp.open(stepOutputFileName, "r")

    for line1 in stepOutFile1:
      line2 = stepOutFile2.readline()
      outFile.write(line1+" "+line2)

    outFile.close()
    stepOutFile1.close()
    stepOutFile2.close()

if os.access("mergerTree.log", os.F_OK): os.unlink("mergerTree.log")
