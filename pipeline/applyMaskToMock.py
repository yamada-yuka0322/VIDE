#!/usr/bin/env python

# applies a mask to a given dataset

import numpy as np
import os
import sys
import void_python_tools as vp
import argparse
import imp
import healpy as hp

# ------------------------------------------------------------------------------

def my_import(name):
    mod = __import__(name)
    components = name.split('.')
    for comp in components[1:]:
        mod = getattr(mod, comp)
    return mod

# -----------------------------------------------------------------------------

LIGHT_SPEED = 299792.458

parser = argparse.ArgumentParser(description='options')
parser.add_argument('--scripts', dest='script', action='store_const',
                   const=True, default=False,
                   help='write scripts')
parser.add_argument('--parmFile', dest='parmFile',
                   default="",
                   help='path to parameter file')
parser.add_argument('--subsamples', dest='subsample', action='store_const',
                   const=True, default=False,
                   help='write subsamples')
parser.add_argument('--halos', dest='halos', action='store_const',
                   const=True, default=False,
                   help='write halos')
parser.add_argument('--hod', dest='hod', action='store_const',
                   const=True, default=False,
                   help='write hod')
parser.add_argument('--all', dest='all', action='store_const',
                   const=True, default=False,
                   help='write everything')
args = parser.parse_args()


filename = args.parmFile
print " Loading parameters from", filename
if not os.access(filename, os.F_OK):
  print "  Cannot find parameter file %s!" % filename
  exit(-1)
parms = imp.load_source("name", filename)
globals().update(vars(parms))


#------------------------------------------------------------------------------
def getSampleName(setName, redshift, useVel, iSlice=-1, iVol=-1):

  sampleName = setName

  sampleName += "_z" + redshift

  if iVol != -1: sampleName += "_d" + iVol

  return sampleName

#------------------------------------------------------------------------------
# for given dataset parameters, outputs a script for use with analyzeVoids
def writeObservationScript(setName, dataFileNameBase, maskFileName,
                scriptDir, catalogDir, fileNums, redshifts, 
                useVel, minRadius, omegaM,
                suffix=".dat"):


  if useVel: setName += "_pv"

  scriptFileName = scriptDir + "/" + setName + ".py"
  scriptFile = open(scriptFileName, 'w')

  scriptFile.write("""#!/usr/bin/env/python
import os
from void_python_tools.backend.classes import *

continueRun = True # set to True to enable restarting aborted jobs
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

dataPortions = ["central", "all"]
dataSampleList = []
""")


  dataInfo = """
setName = "{setName}"

workDir = "{voidOutputDir}/{setName}/"
inputDataDir = "{inputDataDir}"
figDir = "{figDir}/{setName}/"
logDir = "{logDir}/{setName}/"

numZobovDivisions = {numZobovDivisions}
numZobovThreads = {numZobovThreads}
               """
  scriptFile.write(dataInfo.format(setName=setName, figDir=figDir,
                                   logDir=logDir, voidOutputDir=voidOutputDir,
                                   inputDataDir=catalogDir, 
                                   numZobovDivisions=numZobovDivisions,
                                   numZobovThreads=numZobovThreads))

  sampleInfo = """
newSample = Sample(dataFile = "{dataFile}",
                   fullName = "{sampleName}",
                   nickName = "{sampleName}",
                   dataType = "observation",
                   maskFile = "{maskFileName}",
                   zBoundary = ({zMin}, {zMax}),
                   zRange    = ({zMin}, {zMax}),
                   minVoidRadius = {minRadius},
                   includeInHubble = True,
                   partOfCombo = False,
                   isCombo = False,
                   usePecVel = {usePecVel})
dataSampleList.append(newSample)
newSample.addStack({zMin}, {zMax}, {minRadius}  , {minRadius}+2, True, False)
newSample.addStack({zMin}, {zMax}, {minRadius}  , {minRadius}+4, True, False)
newSample.addStack({zMin}, {zMax}, {minRadius}+2, {minRadius}+6, True, False)
newSample.addStack({zMin}, {zMax}, {minRadius}+6, {minRadius}+10, True, False)
newSample.addStack({zMin}, {zMax}, {minRadius}+10, {minRadius}+18, True, False)
newSample.addStack({zMin}, {zMax}, {minRadius}+18, {minRadius}+24, True, False)
               """
  for (iFile, redshift) in enumerate(redshifts):
    fileNum = fileNums[iFile]
    zBox = float(redshift)
   
    zMin = zBox
    zMax = zBox + 0.1 
    zMin, zMax = getMaskedZRange(lbox, zBox, zMin, zMax, omegaM)
    dataFileName = dataFileNameBase + fileNum + suffix

    sampleName = getSampleName(setName, redshift, useVel,
                                iSlice=-1, iVol=-1)

    scriptFile.write(sampleInfo.format(dataFile=dataFileName,
                                         dataFormat=dataFormat,
                                         sampleName=sampleName,
                                         maskFileName=maskFileName,
                                         zMin=zMin,
                                         zMax=zMax,
                                         usePecVel=useVel,
                                         minRadius=minRadius))

  scriptFile.close()
  return


#------------------------------------------------------------------------------
# now place these particles on a lightcone, restrict redshift range, apply mask
def applyMask(inFileName, outFileName, maskFileName, lbox, zBox, zMin, zMax, omegaM):
  mask = hp.read_map(maskFileName)
  nside = hp.get_nside(mask)

  inFile  = open(inFileName, 'r') 
  outFile = open(outFileName, 'w')

  Om = float(omegaM)

  # converter from redshift to comoving distance   
  zVsDY = np.linspace(0., zBox+8*lbox*100./LIGHT_SPEED, 10000)
  zVsDX = np.zeros(len(zVsDY))
  for i in xrange(len(zVsDY)):
    zVsDX[i] = vp.angularDiameter(zVsDY[i], Om=Om)

  zBoxInMpc = vp.angularDiameter(zBox, Om=Om)

  for (iLine,line) in enumerate(inFile):
    if iLine < 5:
      print >> outFile, line.rstrip()
      continue

    line = line.split(' ')
    uniqueID = int(line[0])
    x  = float(line[1]) - lbox/2.
    y  = float(line[2]) - lbox/2.
    z  = float(line[3]) - lbox/2.
    vz = float(line[4])

    # TODO if usePecvel

    redshift = np.sqrt(x*x + y*y + z*z) 
    redshift = np.interp(zBoxInMpc+100./LIGHT_SPEED*redshift, zVsDX, zVsDY)

    if redshift < zMin or redshift > zMax: continue

    vec = np.array((x,y,z))
    theta, phi = hp.vec2ang(vec)
    theta = theta[0]
    phi = phi[0]
    RA = phi*180./np.pi
    Dec = 180./np.pi*(np.pi/2.-theta)

    pix = hp.vec2pix(nside, x, y, z)
    if mask[pix] <= 0.2: continue   
  
    print >> outFile, RA, Dec, redshift*LIGHT_SPEED, uniqueID, x, y, z
 
  inFile.close()
  outFile.close()

#------------------------------------------------------------------------------
def getMaskedZRange(lbox, zBox, zMin, zMax, omegaM):

  if zMin < zBox: zMin = zBox 
  Om = float(omegaM)

  # converter from redshift to comoving distance   
  zVsDY = np.linspace(0., zBox+8*lbox*100./LIGHT_SPEED, 10000)
  zVsDX = np.zeros(len(zVsDY))
  for i in xrange(len(zVsDY)):
    zVsDX[i] = vp.angularDiameter(zVsDY[i], Om=Om)

  zBoxInMpc = vp.angularDiameter(zBox, Om=Om)

  boxWidthZ = np.interp(vp.angularDiameter(zBox,Om=Om)+100. / \
                  LIGHT_SPEED*lbox/2., zVsDX, zVsDY)-zBox

  print "RANGE", zBox, zBox+boxWidthZ
  if zMax > zBox+boxWidthZ: zMax = zBox+boxWidthZ
   
  return zMin, zMax

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

inCatalogDir = catalogDir

prefix = "masked_" + prefix
voidOutputDir += "masked/"
figDir += "masked/"
logDir += "masked/"
scriptDir += "masked/"
catalogDir += "masked/"
dataType = "observation"

if not os.access(scriptDir, os.F_OK): os.mkdir(scriptDir)
if not os.access(catalogDir, os.F_OK): os.mkdir(catalogDir)

#------------------------------------------------------------------------------
# first the directly downsampled runs
# Note: ss0.002   ~ SDSS DR7 dim2
#       ss0.0004 ~ SDSS DR9 mid 
baseResolution = float(numPart)/lbox/lbox/lbox # particles/Mpc^3
for thisSubSample in subSamples:

  keepFraction = float(thisSubSample) / baseResolution
  maxKeep = keepFraction * numPart
  minRadius = int(np.ceil(lbox/maxKeep**(1./3)))

  if args.script or args.all:
    print " Doing subsample", thisSubSample, " scripts"
    setName = prefix+"ss"+str(thisSubSample)
    writeObservationScript(setName, "masked_md.ss"+str(thisSubSample)+"_z", 
                  maskFileName, 
                  scriptDir, catalogDir, fileNums, 
                  redshifts, False, minRadius, omegaM)
    writeObservationScript(setName, "masked_md.ss"+str(thisSubSample)+"_z", 
                  maskFileName, 
                  scriptDir, catalogDir, fileNums, 
                  redshifts, True, minRadius, omegaM)
  
  if args.subsample or args.all:
    print " Doing subsample", thisSubSample

    for (iRedshift, redshift) in enumerate(redshifts):
      print "   redshift", redshift
  
      zMin = float(redshift)
      zMax = zMin + 0.1 
      zMin, zMax = getMaskedZRange(lbox, float(redshift), zMin, zMax, omegaM)

      if dataFormat == "multidark" or dataFormat == "random":
        inFileName = inCatalogDir+"/md.ss"+str(thisSubSample)+"_z"+redshift+".dat"
        outFileName = catalogDir+"/masked_md.ss"+str(thisSubSample)+"_z"+redshift+".dat"
        
        applyMask(inFileName, outFileName, maskFileName, lbox, float(redshift), 
                  zMin, zMax, omegaM)


# -----------------------------------------------------------------------------
# now halos
if (args.script or args.all) and dataFormat == "multidark":
  print " Doing halo scripts"

  for minHaloMass in minHaloMasses:
  
    if dataFormat == "multidark":
      setName = prefix+"halos_min"+str(minHaloMass)
      writeScript(setName, "md.halos_min"+str(minHaloMass)+"_z", 
                scriptDir, catalogDir, fileNums, 
                redshifts, 
                numSubvolumes, numSlices, False, lbox, minRadius, omegaM)
      writeScript(setName, "md.halos_min"+str(minHaloMass)+"_z", 
                scriptDir, catalogDir, fileNums, 
                redshifts, 
                numSubvolumes, numSlices, True, lbox, minRadius, omegaM)

if args.halos or args.all: 
  print " Doing halos"

  for minHaloMass in minHaloMasses:
    print "  min halo mass = ", minHaloMass

    for (iRedshift, redshift) in enumerate(redshifts):
      print "   z = ", redshift

      dataFile = catalogDir+haloFileBase+fileNums[iRedshift]
      inFile = open(dataFile, 'r')
      numPart = 0
      for line in inFile: 
        line = line.split(',')
        if minHaloMass == "none" or float(line[6]) > minHaloMass:
          numPart += 1
      inFile.close()

      sampleName = "md.halos_min"+str(minHaloMass)+"_z"+redshift
      outFile = open(catalogDir+"/"+sampleName+".dat", 'w')

# -----------------------------------------------------------------------------
# now the SDSS HOD

if (args.script or args.all) and dataFormat == "multidark":
  print " Doing DR7 HOD scripts"
  if dataFormat == "multidark":
    setName = prefix+"hod_dr72dim2"
    writeScript(setName, "md.hod_dr72dim2_z",
              scriptDir, catalogDir, fileNums, redshifts, 
              numSubvolumes, numSlices, False, lbox, 5, omegaM)
    writeScript(setName, "md.hod_dr72dim2_z",
              scriptDir, catalogDir, fileNums, redshifts, 
              numSubvolumes, numSlices, True, lbox, 5, omegaM)

if args.hod or args.all:
  print " Doing DR7 HOD"
  for (iRedshift, redshift) in enumerate(redshifts):
    print "  z = ", redshift

    sampleName = getSampleName("md.hod_dr72dim2", redshift, False)
    outFileName = catalogDir+"/"+sampleName+".dat"

# -----------------------------------------------------------------------------
# now the BOSS HOD
if (args.script or args.all) and dataFormat == "multidark":
  print " Doing DR9 HOD scripts"
  if dataFormat == "multidark":
    setName = prefix+"hod_dr9mid"
    writeScript(setName, "md.hod_dr9mid_z",
              scriptDir, catalogDir, fileNums, redshifts, 
               numSubvolumes, numSlices, False, lbox, 15, omegaM)
    writeScript(setName, "md.hod_dr9mid_z",
              scriptDir, catalogDir, fileNums, redshifts, 
               numSubvolumes, numSlices, True, lbox, 15, omegaM)

if args.hod or args.all:
  print " Doing DR9 HOD"
  for (iRedshift, redshift) in enumerate(redshifts):
    print "  z = ", redshift

    outFileName = catalogDir+"/"+sampleName+".dat"

    if dataFormat == "multidark" or dataFormat == "random":
      sampleName = getSampleName("md.hod_dr9mid", redshift, False)
      inFileName = inCatalogDir+"/"+sampleName+".dat"
      outFileName = catalogDir+"/masked_"+sampleName+".dat"
      
      zMin = float(redshift)
      zMax = zMin + 0.1 
      zMin, zMax = getMaskedZRange(lbox, float(redshift), zMin, zMax, omegaM)
        
      applyMask(inFileName, outFileName, maskFileName, lbox, float(redshift), 
                zMin, zMax, omegaM)
