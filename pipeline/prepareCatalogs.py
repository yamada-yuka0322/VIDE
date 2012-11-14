#!/usr/bin/env python

# prepares input catalogs based on multidark simulations
#   (borrows heavily from generateMock, but doesn't hold much in memory)
# also creates necessary analyzeVoids input files

import numpy as np
import os
import sys
import void_python_tools as vp
import argparse

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# CONFIGURATION

# directory for the input simulation/observational particle files
catalogDir = os.getenv("HOME")+"/workspace/Voids/catalogs/multidark/"

# path to HOD code
hodPath = os.getenv("HOME")+"/projects/Voids/hod/HOD.x"

# where to put the final void catalog, figures, and output logs
voidOutputDir = os.getenv("HOME")+"/workspace/Voids/multidark/"
figDir = os.getenv("PWD")+"/../figs/multidark/"
logDir = os.getenv("PWD")+"/../logs/multidark/"

# where to place the pipeline scripts
scriptDir = os.getenv("PWD")+"/multidark/"

# simulation or observation?
dataType = "simulation"

# available formats for simulation: gadget, multidark
dataFormat = "multidark"
dataUnit = 1 # as multiple of Mpc/h

# place particles on the lightcone?
useLightCone = True

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
prefix = "md_"

# list of desired subsamples
subSamples = ((0.1, 0.05, 0.01, 0.002, 0.001, 0.0004, 0.0002))

# simulation information
numPart = 1024*1024*1024
lbox = 1000 # Mpc/h
omegaM = 0.27
hubble = 0.70

# END CONFIGURATION
# -----------------------------------------------------------------------------
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
                   help='write hos')
parser.add_argument('--all', dest='all', action='store_const',
                   const=True, default=False,
                   help='write everything')
args = parser.parse_args()

#------------------------------------------------------------------------------
def getSampleName(setName, redshift, useVel, iSlice=-1, iVol=-1):

  sampleName = setName

  #if useVel: sampleName += "_pv"

  sampleName += "_z" + redshift

  #if iSlice != -1: sampleName += "_s" + str(iSlice)

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
                   dataUnit = {dataUnit},
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
    elif dataFormat == "gadget":
      writeScript(setName, particleFileBase, scriptDir, catalogDir, fileNums,
                   redshifts, 
                  numSubvolumes, numSlices, False, lbox, minRadius, omegaM,
                  subsample=thisSubSample, suffix="")
      writeScript(setName, particleFileBase, scriptDir, catalogDir, fileNums,
                   redshifts, 
                  numSubvolumes, numSlices, True, lbox, minRadius, omegaM,
                  subsample=thisSubSample, suffix="")
  
  if args.subsample or args.all:
    print " Doing subsample", thisSubSample

    for (iRedshift, redshift) in enumerate(redshifts):
      print "   redshift", redshift
      dataFile = catalogDir+"/"+particleFileBase+fileNums[iRedshift]
      inFile = open(dataFile, 'r')

      sampleName = "md.ss"+str(thisSubSample)+"_z"+redshift
      outFile = open(catalogDir+"/"+sampleName+".dat", 'w')

      outFile.write("%f\n" %(lbox))
      outFile.write("%s" %(omegaM))
      outFile.write("%s" %(hubble))
      outFile.write("%s\n" %(redshift))
      outFile.write("%d\n" %(maxKeep))

      numKept = 0
      for line in inFile:
        if np.random.uniform() > keepFraction: continue
        numKept += 1
        if numKept > maxKeep: break
        line = line.split(',')
        x  = float(line[0])
        y  = float(line[1])
        z  = float(line[2])
        vz = float(line[3])

        outFile.write("%e %e %e %e\n" %(x,y,z,vz))

      outFile.write("-99 -99 -99 -99\n")
      print "KEEPING:", numKept, "...predicted:", maxKeep
      inFile.close()
      outFile.close()

# -----------------------------------------------------------------------------
# now halos
if args.script or args.all:
  print " Doing halo scripts"

  dataFile = catalogDir+"/mdr1_halos_z"+fileNums[0]
  inFile = open(dataFile, 'r')
  numPart = 0
  for line in inFile: numPart += 1
  inFile.close()

  minRadius = 2*int(np.ceil(lbox/numPart**(1./3.)))

  if dataFormat == "multidark":
    setName = prefix+"halos"
    writeScript(setName, "md.halos_z", scriptDir, catalogDir, fileNums, 
              redshifts, 
              numSubvolumes, numSlices, False, lbox, minRadius, omegaM)
    writeScript(setName, "md.halos_z", scriptDir, catalogDir, fileNums, 
                redshifts, numSubvolumes, 
                numSlices, True, lbox, minRadius, omegaM)

if args.halos or args.all: 
  print " Doing halos"

  for (iRedshift, redshift) in enumerate(redshifts):
    print "  z = ", redshift

    dataFile = catalogDir+"/mdr1_halos_z"+fileNums[iRedshift]
    inFile = open(dataFile, 'r')
    numPart = 0
    for line in inFile: numPart += 1
    inFile.close()

    sampleName = "md.halos_z"+redshift
    outFile = open(catalogDir+"/"+sampleName+".dat", 'w')

    outFile.write("%f\n" %(lbox))
    outFile.write("%s" %(omegaM))
    outFile.write("%s" %(hubble))
    outFile.write("%s\n" %(redshift))
    outFile.write("%d\n" %(numPart))

    inFile = open(dataFile, 'r')
    numKept = 0
    for line in inFile:
      numKept += 1
      line = line.split(',')
      x  = float(line[0])
      y  = float(line[1])
      z  = float(line[2])
      vz = float(line[5])

      # write to output file
      outFile.write("%e %e %e %e\n" %(x,y,z,vz))

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

if args.script or args.all:
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

    parFileName = "./hod.par"
    parFile = open(parFileName, 'w')
    haloFile = catalogDir+"/mdr1_halos_z"+fileNums[iRedshift]
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

    sampleName = getSampleName("md.hod_dr72dim2", redshift, False)
    outFileName = catalogDir+"/"+sampleName+".dat"
    os.system("mv hod.mock" + " " + outFileName)

    os.system("rm ./hod.*")

# -----------------------------------------------------------------------------
# now the BOSS HOD
if args.script or args.all:
  print " Doing DR9 HOD scripts"
  if dataFormat == "multidark":
    setName = prefix+"hod_dr9mid"
    writeScript(setName, "md.hod_dr9mid_z",
              scriptDir, catalogDir, fileNums, redshifts, 
               numSubvolumes, numSlices, False, lbox, 5, omegaM)
    writeScript(prefix, "md.hod_dr9mid_z",
              scriptDir, catalogDir, fileNums, redshifts, 
               numSubvolumes, numSlices, True, lbox, 5, omegaM)

if args.hod or args.all:
  print " Doing DR9 HOD"
  for (iRedshift, redshift) in enumerate(redshifts):
    print "  z = ", redshift

    parFileName = "./hod.par"
    parFile = open(parFileName, 'w')
    haloFile = catalogDir+"/mdr1_halos_z"+fileNums[iRedshift]
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

    sampleName = getSampleName("md.hod_dr9mid", redshift, False)
    outFileName = catalogDir+"/"+sampleName+".dat"
    os.system("mv hod.mock" + " " + outFileName)

    os.system("rm ./hod.*")


