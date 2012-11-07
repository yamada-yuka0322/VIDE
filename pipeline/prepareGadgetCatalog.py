#!/usr/bin/env python

# script which prepares inputs for gadget-based void run

import numpy as np
import os
import sys
import void_python_tools as vp
import argparse

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# CONFIGURATION

# directory for the input simulation/observational particle files
catalogDir = os.getenv("HOME")+"/workspace/Voids/catalogs/gadget/"

# where to put the final void catalog, figures, and output logs
voidOutputDir = os.getenv("HOME")+"/workspace/Voids/gadget/"
figDir = os.getenv("PWD")+"/../figs/gadget/"
logDir = os.getenv("PWD")+"/../logs/gadget/"

# where to place the pipeline scripts
scriptDir = os.getenv("PWD")+"/gadget/"

# simulation or observation?
dataType = "simulation"

# available formats for simulation: gadget, multidark
dataFormat = "gadget"

# place particles on the lightcone?
useLightCone = False 

# common filename of particle files
particleFileBase = "mdr1_particles_z"

# list of file numbers for the particle files
# to get particle file name, we take particleFileBase+fileNum
fileNums = (("0.0", "0.53", "1.0")) 

# redshift of each file in the above list
redshifts = (("0.0", "0.53", "1.0")) 

# how many independent slices along the z-axis?
numSlices = 4 

# how many subdivisions along the x- and y- axis?
#   ( = 2 will make 4 subvolumes for each slice, = 3 will make 9, etc.)
numSubvolumes = 1 

# prefix to give all outputs
prefix = "gadget_"

# list of desired subsamples
subSamples = [ 1.0, 0.1 ]

# simulation information
numPart = 1024*1024*1024 
lbox = 1000 # Mpc/h
omegaM = 0.27
hubble = 0.70

# END CONFIGURATION
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

#------------------------------------------------------------------------------
LIGHT_SPEED = 299792.458

def getSampleName(setName, redshift, useVel, iSlice=-1, iVol=-1):

  sampleName = setName

  if useVel: setName += "_pv"
    
  if iSlice != -1: sampleName += "_s" + str(iSlice)

  if iVol != -1: sampleName += "_d" + iVol

  return sampleName

#------------------------------------------------------------------------------
# for given dataset parameters, outputs a script for use with analyzeVoids
def writeScript(setName, dataFileNameBase, 
                scriptDir, catalogDir, fileNums, redshifts, numSubvolumes,
                numSlices, useVel, lbox, minRadius, omegaM, subsample=1.0, 
                suffix=".dat"):

  if useVel: setName += "_pv"

  scriptFileName = scriptDir + "/" + setName + ".py"
  scriptFile = open(scriptFileName, 'w')

  scriptFile.write("""#!/usr/bin/env/python
import os
from void_python_tools.backend.classes import *

continueRun = False
startCatalogStage = 1
endCatalogStage   = 4
               
startAPStage = 1
endAPStage = 7

ZOBOV_PATH = os.getenv("PWD")+"/../zobov/"
CTOOLS_PATH = os.getenv("PWD")+"/../c_tools/"
freshStack = True
errorBars = "CALCULATED"
numIncoherentRuns = 100
ranSeed = 101010
useLCDM = False
bias = 1.16

dataPortions = ["central"]
dataSampleList = []
""")

  dataInfo = """
setName = "{setName}"

workDir = "{voidOutputDir}/{setName}/"
inputDataDir = "{inputDataDir}"
figDir = "{figDir}/{setName}/"
logDir = "{logDir}/{setName}/"
               """
  scriptFile.write(dataInfo.format(setName=setName, figDir=figDir, 
                                   logDir=logDir, voidOutputDir=voidOutputDir,
                                   inputDataDir=catalogDir))

  sampleInfo = """
newSample = Sample(dataFile = "{dataFile}",
                   dataFormat = "{dataFormat}",
                   fullName = "{sampleName}",
                   nickName = "{sampleName}",
                   dataType = "simulation",
                   zBoundary = ({zMin}, {zMax}),
                   zRange    = ({zMin}, {zMax}),
                   zBoundaryMpc = ({zMinMpc}, {zMaxMpc}),
                   omegaM    = {omegaM},
                   minVoidRadius = {minRadius},
                   includeInHubble = True,
                   partOfCombo = False,
                   isCombo = False,
                   boxLen = {boxLen},
                   usePecVel = {usePecVel},
                   numSubvolumes = {numSubvolumes},
                   mySubvolume = "{mySubvolume}",
                   numSubDivisions = 4,
                   useLightCone = {useLightCone},
                   subsample = {subsample})
dataSampleList.append(newSample)
               """

  for (iFile, redshift) in enumerate(redshifts):
    fileNum = fileNums[iFile]
    zBox = float(redshift)
    Om = float(omegaM)
    zBoxMpc = LIGHT_SPEED/100.*vp.angularDiameter(zBox, Om=Om)
    boxMaxMpc = zBoxMpc + lbox

    # converter from redshift to comoving distance   
    zVsDY = np.linspace(0., zBox+8*lbox*100./LIGHT_SPEED, 10000)  
    zVsDX = np.zeros(len(zVsDY))
    for i in xrange(len(zVsDY)):
      zVsDX[i] = vp.angularDiameter(zVsDY[i], Om=Om)

    if useLightCone:
      boxWidthZ = np.interp(vp.angularDiameter(zBox,Om=Om)+100. / \
                  LIGHT_SPEED*lbox, zVsDX, zVsDY)-zBox
      dzSafe = 0.03
    else:
      boxWidthZ = lbox*100./LIGHT_SPEED
      dzSafe = 0.0

    for iSlice in xrange(numSlices):
      sliceMin = zBox + dzSafe + iSlice*(boxWidthZ-dzSafe)/numSlices
      sliceMax = zBox + dzSafe + (iSlice+1)*(boxWidthZ-dzSafe)/numSlices

      sliceMinMpc = sliceMin*LIGHT_SPEED/100.
      sliceMaxMpc = sliceMax*LIGHT_SPEED/100.
  
      sliceMin = "%0.2f" % sliceMin
      sliceMax = "%0.2f" % sliceMax
      sliceMinMpc = "%0.1f" % sliceMinMpc
      sliceMaxMpc = "%0.1f" % sliceMaxMpc

      dataFileName = dataFileNameBase + fileNum + suffix

      for iX in xrange(numSubvolumes):
        for iY in xrange(numSubvolumes):
  
          mySubvolume = "%d%d" % (iX, iY)

          sampleName = getSampleName(setName, redshift, useVel, 
                                     iSlice=iSlice, iVol=mySubvolume)

          scriptFile.write(sampleInfo.format(dataFile=dataFileName, 
                                         dataFormat=dataFormat,
                                         sampleName=sampleName,
                                         zMin=sliceMin,
                                         zMax=sliceMax,
                                         zMinMpc=sliceMinMpc,
                                         zMaxMpc=sliceMaxMpc,
                                         omegaM=Om,
                                         boxLen=lbox,
                                         usePecVel=useVel,
                                         minRadius=minRadius,
                                         numSubvolumes=numSubvolumes,
                                         mySubvolume=mySubvolume,
                                         useLightCone=useLightCone,
                                         subsample=subsample))

  scriptFile.close()
  return


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
if not os.access(scriptDir, os.F_OK): os.mkdir(scriptDir)

#------------------------------------------------------------------------------
# first the directly downsampled runs
# Note: ss0.002   ~ SDSS DR7 dim2
#       ss0.0004 ~ SDSS DR9 mid 
baseResolution = float(numPart)/lbox/lbox/lbox # particles/Mpc^3
for thisSubSample in subSamples:

  keepFraction = float(thisSubSample) / baseResolution
  maxKeep = keepFraction * numPart
  minRadius = int(np.ceil(lbox/maxKeep**(1./3)))

  print " Doing subsample", thisSubSample, " scripts"
  setName = prefix+"ss"+str(thisSubSample)
  writeScript(setName, particleFileBase, scriptDir, catalogDir, fileNums, 
                 redshifts, 
                 numSubvolumes, numSlices, True, lbox, minRadius, omegaM,
                 subsample=thisSubSample, suffix="")
  writeScript(setName, particleFileBase, scriptDir, catalogDir, fileNums, 
                 redshifts, 
                 numSubvolumes, numSlices, False, lbox, minRadius, omegaM,
                 subsample=thisSubSample, suffix="")
  
