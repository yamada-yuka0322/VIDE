#!/usr/bin/env python

# does full void analysis. Also generates 2d/1d stacked plots and hubble diagram

from void_python_tools.backend import *
from void_python_tools.plotting import *
import imp
import pickle

# ------------------------------------------------------------------------------

if (len(sys.argv) == 1):
  print "Usage: ./analyzeVoids.py parameter_file.py"
  exit(-1)

if (len(sys.argv) > 1):
  filename = sys.argv[1]
  print " Loading parameters from", filename
  if not os.access(filename, os.F_OK):
    print "  Cannot find parameter file %s!" % filename
    exit(-1)
  parms = imp.load_source("name", filename)
  regenerateFlag = False
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

  fp = open(zobovDir+"/sample_info.txt", 'w')
  fp.write("Sample name: %s\n" % sample.fullName)
  fp.write("Sample nickname: %s\n" % sample.nickName)
  fp.write("Data type: %s\n" % sample.dataType)
  fp.write("Redshift range: %f - %f\n" %(sample.zBoundary[0],sample.zBoundary[1]))
  fp.write("Estimated mean particle separation: %g\n" % sample.minVoidRadius)

  if (sample.dataType == "simulation"):
    fp.write("Particles placed on lightcone: %g\n" % sample.useLightCone)
    fp.write("Peculiar velocities included: %g\n" % sample.usePecVel)
    fp.write("Additional subsampling fraction: %s\n" % sample.subsample[-1])
    fp.write("Simulation box length (Mpc/h): %g\n" % sample.boxLen)
    fp.write("Simulation Omega_M: %g\n" % sample.omegaM)
    fp.write("Number of simulation subvolumes: %s\n" % sample.numSubvolumes)
    fp.write("My subvolume index: %s\n" % sample.mySubvolume)
  fp.close()

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
                   continueRun=continueRun, regenerate=regenerateFlag)

  # --------------------------------------------------------------------------
  if (startCatalogStage <= 2) and (endCatalogStage >= 2) and not sample.isCombo:
    print "  Extracting voids with ZOBOV...",
    sys.stdout.flush()

    launchZobov(sample, ZOBOV_PATH, zobovDir=zobovDir, logDir=logDir, 
                continueRun=continueRun, numZobovDivisions=numZobovDivisions,
                 numZobovThreads=numZobovThreads)

  # -------------------------------------------------------------------------
  if (startCatalogStage <= 3) and (endCatalogStage >= 3) and not sample.isCombo:

    print "  Taking data portions", "...",
    sys.stdout.flush()
        
    logFile = logDir+"/pruneVoids_"+sampleName+".out"

    PRUNE_PATH = CTOOLS_PATH+"/stacking/pruneVoids"

    launchPrune(sample, PRUNE_PATH, 
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
