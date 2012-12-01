#!/usr/bin/env python

# prepares input catalogs based on multidark simulations
#   (borrows heavily from generateMock, but doesn't hold much in memory)
# also creates necessary analyzeVoids input files

import numpy as np
import os
import sys
import void_python_tools as vp
import argparse
import imp

# -----------------------------------------------------------------------------

LIGHT_SPEED = 299792.458

parser = argparse.ArgumentParser(description='options')
parser.add_argument('--scripts', dest='script', action='store_const',
                   const=True, default=False,
                   help='write scripts')
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
parser.add_argument('--parmFile', dest='parmFile',
                   default="",
                   help='path to parameter file')
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

continueRun = False # set to True to enable restarting aborted jobs
startCatalogStage = 1
endCatalogStage   = 4
               
startAPStage = 1
endAPStage = 7

ZOBOV_PATH = "@CMAKE_BINARY_DIR@/zobov/"
CTOOLS_PATH = "@CMAKE_BINARY_DIR@/c_tools/"
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
                   dataFormat = "{dataFormat}",
                   dataUnit = {dataUnit},
                   fullName = "{sampleName}",
                   nickName = "{sampleName}",
                   dataType = "simulation",
                   zBoundary = ({zMin}, {zMax}),
                   zRange    = ({zMin}, {zMax}),
                   zBoundaryMpc = ({zMinMpc}, {zMaxMpc}),
                   omegaM    = {omegaM},
                   minVoidRadius = {minRadius},
                   profileBinSize = 1.0,
                   includeInHubble = True,
                   partOfCombo = False,
                   isCombo = False,
                   boxLen = {boxLen},
                   usePecVel = {usePecVel},
                   numSubvolumes = {numSubvolumes},
                   mySubvolume = "{mySubvolume}",
                   useLightCone = {useLightCone},
                   subsample = {subsample})
dataSampleList.append(newSample)
newSample.addStack({zMin}, {zMax}, 2*{minRadius}  , 2*{minRadius}+2, True, False)
newSample.addStack({zMin}, {zMax}, 2*{minRadius}  , 2*{minRadius}+4, True, False)
newSample.addStack({zMin}, {zMax}, 2*{minRadius}+2, 2*{minRadius}+6, True, False)
newSample.addStack({zMin}, {zMax}, 2*{minRadius}+6, 2*{minRadius}+10, True, False)
newSample.addStack({zMin}, {zMax}, 2*{minRadius}+10, 2*{minRadius}+18, True, False)
newSample.addStack({zMin}, {zMax}, 2*{minRadius}+18, 2*{minRadius}+24, True, False)
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

          sampleName = getSampleName(setName, sliceMin, useVel,
                                     iSlice=iSlice, iVol=mySubvolume)

          scriptFile.write(sampleInfo.format(dataFile=dataFileName,
                                         dataFormat=dataFormat,
                                         dataUnit=dataUnit,
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
    if dataFormat == "multidark":
      writeScript(setName, "md.ss"+str(thisSubSample)+"_z", 
                  scriptDir, catalogDir, fileNums, 
                  redshifts, 
                  numSubvolumes, numSlices, False, lbox, minRadius, omegaM, 
                  subsample=1.0)
      writeScript(setName, "md.ss"+str(thisSubSample)+"_z", 
                  scriptDir, catalogDir, 
                  fileNums, 
                  redshifts, numSubvolumes, numSlices, True, lbox, minRadius, 
                  omegaM, subsample=1.0)
    elif dataFormat == "gadget" or dataFormat == "lanl":
      writeScript(setName, particleFileBase, scriptDir, catalogDir, fileNums,
                   redshifts, 
                  numSubvolumes, numSlices, False, lbox, minRadius, omegaM,
                  subsample=thisSubSample, suffix="")
      writeScript(setName, particleFileBase, scriptDir, catalogDir, fileNums,
                   redshifts, 
                  numSubvolumes, numSlices, True, lbox, minRadius, omegaM,
                  subsample=thisSubSample, suffix="")
    elif dataFormat == "random":
      writeScript(setName, "ran.ss"+str(thisSubSample)+"_z", 
                  scriptDir, catalogDir, fileNums, 
                  redshifts, 
                  numSubvolumes, numSlices, False, lbox, minRadius, omegaM, 
                  subsample=1.0)
  
  if args.subsample or args.all:
    print " Doing subsample", thisSubSample

    for (iRedshift, redshift) in enumerate(redshifts):
      print "   redshift", redshift
  
      if dataFormat == "multidark":
        dataFile = catalogDir+"/"+particleFileBase+fileNums[iRedshift]
        inFile = open(dataFile, 'r')

        sampleName = "md.ss"+str(thisSubSample)+"_z"+redshift
        outFile = open(catalogDir+"/"+sampleName+".dat", 'w')

        outFile.write("%f\n" %(lbox))
        outFile.write("%s\n" %(omegaM))
        outFile.write("%s\n" %(hubble))
        outFile.write("%s\n" %(redshift))
        outFile.write("%d\n" %(maxKeep))

        numKept = 0
        for (i,line) in enumerate(inFile):
          if np.random.uniform() > keepFraction: continue
          numKept += 1
          if numKept > maxKeep: break
          line = line.split(',')
          x  = float(line[0])
          y  = float(line[1])
          z  = float(line[2])
          vz = float(line[3])

          outFile.write("%d %e %e %e %e\n" %(i,x,y,z,vz))

        outFile.write("-99 -99 -99 -99 -99\n")
        inFile.close()
        outFile.close()

      elif dataFormat == "random":
        sampleName = "ran.ss"+str(thisSubSample)+"_z"+redshift
        outFile = open(catalogDir+"/"+sampleName+".dat", 'w')

        outFile.write("%f\n" %(lbox))
        outFile.write("%s\n" %(omegaM))
        outFile.write("%s\n" %(hubble))
        outFile.write("%s\n" %(redshift))
        outFile.write("%d\n" %(maxKeep))

        for i in xrange(int(maxKeep)):
          x  = np.random.uniform()*lbox
          y  = np.random.uniform()*lbox
          z  = np.random.uniform()*lbox

          outFile.write("%d %e %e %e %e\n" % (i, x,y,z, 0.))

        outFile.write("-99 -99 -99 -99 -99\n")
        outFile.close()


# -----------------------------------------------------------------------------
# now halos
if (args.script or args.all) and (dataFormat == "multidark" or dataFormat == "lanl"):
  print " Doing halo scripts"

  for minHaloMass in minHaloMasses:
    # estimate number of halos to get density
    dataFile = catalogDir+haloFileBase+fileNums[0]
    inFile = open(dataFile, 'r')
    numPart = 0
    if dataFormat == "multidark":
      for line in inFile: 
        line = line.split(',')
        if minHaloMass == "none" or float(line[6]) > minHaloMass:
          numPart += 1
    elif dataFormat == "lanl":
      for line in inFile: 
        line = line.split(' ')
        if minHaloMass == "none" or float(line[0]) > minHaloMass:
          numPart += 1
    inFile.close()

    minRadius = 2*int(np.ceil(lbox/numPart**(1./3.)))
  
    setName = prefix+"halos_min"+str(minHaloMass)
    writeScript(setName, prefix+"halos_min"+str(minHaloMass)+"_z", 
                scriptDir, catalogDir, fileNums, 
                redshifts, 
                numSubvolumes, numSlices, False, lbox, minRadius, omegaM)
    writeScript(setName, prefix+"halos_min"+str(minHaloMass)+"_z", 
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
      if dataFormat == "multidark":
        for line in inFile: 
          line = line.split(',')
          if minHaloMass == "none" or float(line[6]) > minHaloMass:
            numPart += 1
      elif dataFormat == "lanl":
          line = line.split(' ')
          if minHaloMass == "none" or float(line[0]) > minHaloMass:
            numPart += 1
      inFile.close()

      sampleName = prefix+"halos_min"+str(minHaloMass)+"_z"+redshift
      outFile = open(catalogDir+"/"+sampleName+".dat", 'w')

      outFile.write("%f\n" %(lbox))
      outFile.write("%s\n" %(omegaM))
      outFile.write("%s\n" %(hubble))
      outFile.write("%s\n" %(redshift))
      outFile.write("%d\n" %(numPart))

      inFile = open(dataFile, 'r')
      for (iHalo,line) in enumerate(inFile):
        if dataFormat == "multidark":
          line = line.split(',')
          if minHaloMass == "none" or float(line[6]) > minHaloMass:
            x  = float(line[0])
            y  = float(line[1])
            z  = float(line[2])
            vz = float(line[5])
        elif dataFormat == "lanl":
          line = line.split(' ')
          if minHaloMass == "none" or float(line[0]) > minHaloMass:
            x  = float(line[1])
            y  = float(line[2])
            z  = float(line[3])
            vz = float(line[4])

          # write to output file
          outFile.write("%d %e %e %e %e\n" %(iHalo,x,y,z,vz))

      outFile.write("-99 -99 -99 -99 -99\n")
      inFile.close()
      outFile.close()

# -----------------------------------------------------------------------------
# now the SDSS HOD
parFileText = """
% cosmology
OMEGA_M {omegaM}
HUBBLE  {hubble}
OMEGA_B 0.0469
SIGMA_8 0.82
SPECTRAL_INDX 0.95
ITRANS 5
REDSHIFT {redshift}

% halo definition
%DELTA_HALO 200
DELTA_HALO 740.74 % 200/Om_m
M_max 1.00E+16

% fit function types
pdfs 11
pdfc 2
EXCLUSION 4

% hod parameters 
M_min {Mmin}
GALAXY_DENSITY 0.0111134 % computed automatically if M_min set, use for sanity
M1 {M1}   
sigma_logM {sigma_logM}
alpha {alpha}
M_cut {Mcut}

% simulation info
real_space_xi 1
HOD 1
populate_sim 1
HaloFile {haloFile}
RESOLUTION {numPartPerSide}
BOX_SIZE   {boxSize}

% output
root_filename hod
               """

if (args.script or args.all) and (dataFormat == "multidark" or dataFormat == "lanl"):
  print " Doing DR7 HOD scripts"
  setName = prefix+"hod_dr72dim2"
  writeScript(setName, prefix+"hod_dr72dim2_z",
              scriptDir, catalogDir, fileNums, redshifts, 
              numSubvolumes, numSlices, False, lbox, 5, omegaM)
  writeScript(setName, prefix+"hod_dr72dim2_z",
              scriptDir, catalogDir, fileNums, redshifts, 
              numSubvolumes, numSlices, True, lbox, 5, omegaM)

if args.hod or args.all:
  print " Doing DR7 HOD"
  for (iRedshift, redshift) in enumerate(redshifts):
    print "  z = ", redshift

    parFileName = "./hod.par"
    parFile = open(parFileName, 'w')
    haloFile = catalogDir+haloFileBase+fileNums[iRedshift]
    parFile.write(parFileText.format(omegaM=omegaM,
                                     hubble=hubble,
                                     redshift=redshift,
                                     Mmin=1.99526e12,
                                     M1=3.80189e13,
                                     sigma_logM=0.21,
                                     alpha=1.12,
                                     Mcut=6.91831e11,
                                     haloFile=haloFile,
                                     numPartPerSide=numPart**(1/3.),
                                     boxSize=lbox))
    parFile.close()

    os.system(hodPath+" "+parFileName+">& /dev/null")

    sampleName = getSampleName(prefix+"hod_dr72dim2", redshift, False)
    outFileName = catalogDir+"/"+sampleName+".dat"
    os.system("mv hod.mock" + " " + outFileName)

    os.system("rm ./hod.*")

# -----------------------------------------------------------------------------
# now the BOSS HOD
if (args.script or args.all) and (dataFormat == "multidark" or dataFormat == "lanl"):
  print " Doing DR9 HOD scripts"
  setName = prefix+"hod_dr9mid"
  writeScript(setName, prefix+"hod_dr9mid_z",
              scriptDir, catalogDir, fileNums, redshifts, 
               numSubvolumes, numSlices, False, lbox, 15, omegaM)
  writeScript(setName, prefix+"hod_dr9mid_z",
              scriptDir, catalogDir, fileNums, redshifts, 
               numSubvolumes, numSlices, True, lbox, 15, omegaM)

if args.hod or args.all:
  print " Doing DR9 HOD"
  for (iRedshift, redshift) in enumerate(redshifts):
    print "  z = ", redshift

    # these parameters come from Manera et al 2012, eq. 26
    parFileName = "./hod.par"
    parFile = open(parFileName, 'w')
    haloFile = catalogDir+haloFileBase+fileNums[iRedshift]
    parFile.write(parFileText.format(omegaM=omegaM,
                                     hubble=hubble,
                                     redshift=redshift,
                                     Mmin=1.23e13,
                                     M1=1.e14,
                                     sigma_logM=0.596,
                                     alpha=1.0127,
                                     Mcut=1.19399e13,
                                     haloFile=haloFile,
                                     numPartPerSide=numPart**(1/3.),
                                     boxSize=lbox))
    parFile.close()

    os.system(hodPath+" "+parFileName+">& /dev/null")

    sampleName = getSampleName(prefix+"hod_dr9mid", redshift, False)
    outFileName = catalogDir+"/"+sampleName+".dat"
    os.system("mv hod.mock" + " " + outFileName)

    os.system("rm ./hod.*")


