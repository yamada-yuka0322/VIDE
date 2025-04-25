#!/usr/bin/env python
#+
#   VIDE -- Void IDentification and Examination -- ./python_tools/void_pipeline/__init__.py
#   Copyright (C) 2010-2014 Guilhem Lavaux
#   Copyright (C) 2011-2014 P. M. Sutter
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

# Stage 1 : generate particles
# Stage 2 : find voids
# Stage 3 : prune catalog

from vide.backend import *
import vide
import imp
import os
import pickle

# ------------------------------------------------------------------------------

if (len(sys.argv) == 1):
  print("Usage: ./generateCatalog.py parameter_file.py")
  exit(-1)

if (len(sys.argv) > 1):
  filename = sys.argv[1]
  print(" Loading parameters from", filename)
  if not os.access(filename, os.F_OK):
    print("  Cannot find parameter file %s!" % filename)
    exit(-1)
  parms = imp.load_source("name", filename)
  regenerateFlag = False
  globals().update(vars(parms))
  void_path = os.path.split(vide.__file__)[0]
  ZOBOV_PATH=f'{void_path}/bin/'
  CTOOLS_PATH=ZOBOV_PATH
  print(f"ZOBOV_PATH is  {ZOBOV_PATH}")
else:
  print(" Using default parameters")

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

if not os.access(logDir, os.F_OK):
  os.makedirs(logDir)

if not os.access(figDir, os.F_OK):
  os.makedirs(figDir)

if not continueRun:
  print(" Cleaning out log files...")

  if startCatalogStage <= 1 and glob.glob(logDir+"/generate*") != []:
    os.system("rm %s/generate*" % logDir)
  if startCatalogStage <= 2 and glob.glob(logDir+"/zobov*") != []:
    os.system("rm %s/zobov*" % logDir)
  if startCatalogStage <= 3 and glob.glob(logDir+"/prune*") != []:
    os.system("rm %s/prune*" % logDir)

for sample in dataSampleList:

  sampleName = sample.fullName

  print(" Working with data set", sampleName, "...")
  zobovDir = workDir+"/sample_"+sampleName+"/"
  sample.zobovDir = zobovDir

  if not os.access(zobovDir, os.F_OK):
    os.makedirs(zobovDir)

# ---------------------------------------------------------------------------
  if (startCatalogStage <= 1) and (endCatalogStage >= 1) and not sample.isCombo:
    print("  Extracting tracers from catalog...", end=' ')
    sys.stdout.flush()

    logFile = logDir+"/generate_"+sampleName+".out"

    if sample.dataType == "observation":
      GENERATE_PATH = CTOOLS_PATH+"/generateFromCatalog"
    else:
      GENERATE_PATH = CTOOLS_PATH+"/generateMock"

    launchGenerate(sample, GENERATE_PATH, workDir=workDir, 
                   inputDataDir=inputDataDir, zobovDir=zobovDir,
                   figDir=figDir, logFile=logFile, useComoving=sample.useComoving,
                   continueRun=continueRun, regenerate=regenerateFlag)

  # --------------------------------------------------------------------------
  if (startCatalogStage <= 2) and (endCatalogStage >= 2) and not sample.isCombo:
    print("  Extracting voids with ZOBOV...", end=' ')
    sys.stdout.flush()

    launchZobov(sample, ZOBOV_PATH, zobovDir=zobovDir, logDir=logDir, 
                continueRun=continueRun, numZobovDivisions=numZobovDivisions,
                 numZobovThreads=numZobovThreads, mergingThreshold=mergingThreshold)

  # -------------------------------------------------------------------------
  if (startCatalogStage <= 3) and (endCatalogStage >= 3) and not sample.isCombo:

    print("  Taking data portions", "...", end=' ')
    sys.stdout.flush()
        
    logFile = logDir+"/pruneVoids_"+sampleName+".out"

    PRUNE_PATH = CTOOLS_PATH+"/pruneVoids"

    launchPrune(sample, PRUNE_PATH, 
                logFile=logFile, zobovDir=zobovDir, 
                useComoving=sample.useComoving, continueRun=continueRun, mergingThreshold=mergingThreshold)

# -------------------------------------------------------------------------
  if (startCatalogStage <= 4) and (endCatalogStage >= 4) and not sample.isCombo:

    print("  Cleaning voids...", end=' ')
    sys.stdout.flush()

    Output(zobovdir=zobovDir, samplename=sampleName)

  #for thisDataPortion in dataPortions:
    #plotRedshiftDistribution(workDir, dataSampleList, figDir, showPlot=False, 
    #                         dataPortion=thisDataPortion, setName=setName)
    #plotSizeDistribution(workDir, dataSampleList, figDir, showPlot=False, 
    #                         dataPortion=thisDataPortion, setName=setName)
    #plotNumberDistribution(workDir, dataSampleList, figDir, showPlot=False, 
    #                         dataPortion=thisDataPortion, setName=setName)
    #plotVoidDistribution(workDir, dataSampleList, figDir, showPlot=False, 
    #                         dataPortion=thisDataPortion, setName=setName)

print("\n Done!")
