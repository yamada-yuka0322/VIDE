#!/usr/bin/env python

# prepares input catalogs based on multidark simulations
#   (borrows heavily from generateMock, but doesn't hold much in memory)
# also creates necessary analyzeVoids input files

import numpy as np
import os
import sys
from Scientific.IO.NetCDF import NetCDFFile
import void_python_tools as vp
import argparse

catalogDir = os.getenv("HOME")+"/workspace/Voids/catalogs/multidark/"
hodPath = os.getenv("HOME")+"/projects/Voids/hod/HOD.x"

voidOutputDir = os.getenv("HOME")+"/workspace/Voids/multidark/"
scriptDir = os.getenv("HOME")+"/projects/Voids/scripts/multidark/"
figDir = os.getenv("HOME")+"/projects/Voids/figs/multidark/"
logDir = os.getenv("HOME")+"/projects/Voids/logs/multidark/"

dataType = "simulation"
dataFormat = "gadget"
useLightCone = False # place particles on a light cone?

redshiftFileBase = "mdr1_particles_z" # common filename of particle files
redshifts = (("0.0", "0.53", "1.0")) # list of redshift particle files

numSlices = 4 # how many slices along the z-axis?
numSubvolumes = 1 # how many subdivisions along the x- and y-axes?

prefix = "testrun_" # prefix to give all outputs
subSamples = [ 0.01 ] # list of desired subsamples

numPart = 1024*1024*1024 # number of particles in base catalog
lbox = 1000 # Mpc/h
omegaM = 0.27
hubble = 0.70

# -----------------------------------------------------------------------------
LIGHT_SPEED = 299792.458

#------------------------------------------------------------------------------
def getSampleName(setName, redshift, useVel, iSlice=-1, iVol=-1):

  sampleName = setName

  if useVel: setName += "_pv"
    
  if iSlice != -1: sampleName += "_s" + str(iSlice)

  if iVol != -1: sampleName += "_d" + iVol

  return sampleName

#------------------------------------------------------------------------------
# for given dataset parameters, outputs a script for use with analyzeVoids
def writeScript(setName, dataFileNameBase, 
                scriptDir, catalogDir, redshifts, numSubvolumes,
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
endCatalogStage   = 3
               
startAPStage = 1
endAPStage = 6

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

  for redshift in redshifts:
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

      dataFileName = dataFileNameBase + redshift + suffix

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
baseResolution = numPart/lbox/lbox/lbox # particles/Mpc^3
for thisSubSample in subSamples:

  keepFraction = float(thisSubSample) / baseResolution
  maxKeep = keepFraction * numPart
  minRadius = int(np.ceil(lbox/maxKeep**(1./3)))

  print " Doing subsample", thisSubSample, " scripts"
  setName = prefix+"ss"+str(thisSubSample)
  writeScript(setName, redshiftFileBase, scriptDir, catalogDir, redshifts, 
                 numSubvolumes, numSlices, True, lbox, minRadius, omegaM,
                 subsample=thisSubSample, suffix="")
  writeScript(setName, redshiftFileBase, scriptDir, catalogDir, redshifts, 
                 numSubvolumes, numSlices, False, lbox, minRadius, omegaM,
                 subsample=thisSubSample, suffix="")
  
