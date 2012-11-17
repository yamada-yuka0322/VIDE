#!/usr/bin/env python

# does full void analysis. Also generates 2d/1d stacked plots and hubble diagram

from void_python_tools.backend import *
from void_python_tools.plotting import *
import imp
import pickle

# ------------------------------------------------------------------------------

def my_import(name):
    mod = __import__(name)
    components = name.split('.')
    for comp in components[1:]:
        mod = getattr(mod, comp)
    return mod

if (len(sys.argv) == 0):
  print "Usage: ./anylyzeVoids.py parameter_file.py"
  exit(-1)

if (len(sys.argv) > 1):
  filename = sys.argv[1]
  print " Loading parameters from", filename
  if not os.access(filename, os.F_OK):
    print "  Cannot find parameter file %s!" % filename
    exit(-1)
  parms = imp.load_source("name", filename)
  globals().update(vars(parms))
else:
  print " Using default parameters"

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

if not os.access(logDir, os.F_OK):
  os.makedirs(logDir)

if not os.access(figDir, os.F_OK):
  os.makedirs(figDir)

if not continueRun:
  print " Cleaning out log files..."

  if startCatalogStage <= 1 and glob.glob(logDir+"/generate*") != []:
    os.system("rm %s/generate*" % logDir)
  if startCatalogStage <= 2 and glob.glob(logDir+"/zobov*") != []:
    os.system("rm %s/zobov*" % logDir)
  if startCatalogStage <= 3 and glob.glob(logDir+"/prune*") != []:
    os.system("rm %s/prune*" % logDir)

for sample in dataSampleList:

  sampleName = sample.fullName

  print " Working with data set", sampleName, "..."
  zobovDir = workDir+"/sample_"+sampleName+"/"
  sample.zobovDir = zobovDir

  if not os.access(zobovDir, os.F_OK):
    os.makedirs(zobovDir)

  # save this sample's information
  with open(zobovDir+"/sample_info.dat", 'w') as output:
    pickle.dump(sample, output, pickle.HIGHEST_PROTOCOL)
  #fp = open(zobovDir+"/sample_info.txt", 'w') 
  #fp.write("Redshift range: %f - %f" %(sample.zBoundary[0], sample.zBoundary[1])
  #fp.close()

# ---------------------------------------------------------------------------
  if (startCatalogStage <= 1) and (endCatalogStage >= 1) and not sample.isCombo:
    print "  Extracting tracers from catalog...",
    sys.stdout.flush()

    logFile = logDir+"/generate_"+sampleName+".out"

    if sample.dataType == "observation":
      GENERATE_PATH = CTOOLS_PATH+"/mock/generateFromCatalog"
    else:
      GENERATE_PATH = CTOOLS_PATH+"/mock/generateMock"

    launchGenerate(sample, GENERATE_PATH, workDir=workDir, 
                   inputDataDir=inputDataDir, zobovDir=zobovDir,
                   figDir=figDir, logFile=logFile, useLCDM=useLCDM,
                   continueRun=continueRun)

  # --------------------------------------------------------------------------
  if (startCatalogStage <= 2) and (endCatalogStage >= 2) and not sample.isCombo:
    print "  Extracting voids with ZOBOV...",
    sys.stdout.flush()

    launchZobov(sample, ZOBOV_PATH, zobovDir=zobovDir, logDir=logDir, 
                continueRun=continueRun, numZobovDivisions=numZobovDivisions,
                 numZobovThreads=numZobovThreads)

  # -------------------------------------------------------------------------
  if (startCatalogStage <= 3) and (endCatalogStage >= 3) and not sample.isCombo:

    for thisDataPortion in dataPortions:
      print "  Taking data portion:", thisDataPortion, "...",
      sys.stdout.flush()
        
      logFile = logDir+"/pruneVoids_"+sampleName+"_"+\
                thisDataPortion+".out"

      PRUNE_PATH = CTOOLS_PATH+"/stacking/pruneVoids"

      launchPrune(sample, PRUNE_PATH, thisDataPortion=thisDataPortion, 
                  logFile=logFile, zobovDir=zobovDir, continueRun=continueRun)

# -------------------------------------------------------------------------
if (startCatalogStage <= 4) and (endCatalogStage >= 4):

  print "  Plotting...",
  sys.stdout.flush()

  for thisDataPortion in dataPortions:
    plotRedshiftDistribution(workDir, dataSampleList, figDir, showPlot=False, 
                             dataPortion=thisDataPortion, setName=setName)
    plotSizeDistribution(workDir, dataSampleList, figDir, showPlot=False, 
                             dataPortion=thisDataPortion, setName=setName)
    plotNumberDistribution(workDir, dataSampleList, figDir, showPlot=False, 
                             dataPortion=thisDataPortion, setName=setName)

print "\n Done!"
