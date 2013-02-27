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
def writeScript(setName, dataFileNameBase, dataFormat,
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
endCatalogStage   = 3
               
startAPStage = 1
endAPStage = 1

regenerateFlag = False
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
                   profileBinSize = "auto",
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
newSample.addStack(0.0, 5.0, 20, 25, False, False)
newSample.addStack(0.0, 5.0, 30, 35, False, False)
newSample.addStack(0.0, 5.0, 40, 45, False, False)
newSample.addStack(0.0, 5.0, 50, 55, False, False)
newSample.addStack(0.0, 5.0, 60, 65, False, False)
newSample.addStack(0.0, 5.0, 70, 75, False, False)
newSample.addStack(0.0, 5.0, 80, 85, False, False)
newSample.addStack(0.0, 5.0, 90, 95, False, False)
#newSample.addStack({zMin}, {zMax}, 2*{minRadius}  , 2*{minRadius}+2, True, False)
#newSample.addStack({zMin}, {zMax}, 2*{minRadius}  , 2*{minRadius}+4, True, False)
#newSample.addStack({zMin}, {zMax}, 2*{minRadius}+2, 2*{minRadius}+6, True, False)
#newSample.addStack({zMin}, {zMax}, 2*{minRadius}+6, 2*{minRadius}+10, True, False)
#newSample.addStack({zMin}, {zMax}, 2*{minRadius}+10, 2*{minRadius}+18, True, False)
#newSample.addStack({zMin}, {zMax}, 2*{minRadius}+18, 2*{minRadius}+24, True, False)
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

    boxWidthZ = np.interp(vp.angularDiameter(zBox,Om=Om)+100. / \
                LIGHT_SPEED*lbox, zVsDX, zVsDY)-zBox

    for iSlice in xrange(numSlices):

      if useLightCone:
        dzSafe = 0.03
        sliceMin = zBox + dzSafe + iSlice*(boxWidthZ-2.*dzSafe)/numSlices
        sliceMax = zBox + dzSafe + (iSlice+1)*(boxWidthZ-2.*dzSafe)/numSlices
        sliceMinMpc = sliceMin*LIGHT_SPEED/100.
        sliceMaxMpc = sliceMax*LIGHT_SPEED/100.
      else:
        sliceMinMpc = zBoxMpc + iSlice*lbox/numSlices
        sliceMaxMpc = zBoxMpc + (iSlice+1)*lbox/numSlices
        sliceMin = np.interp(sliceMinMpc*100./LIGHT_SPEED, zVsDX, zVsDY)
        sliceMax = np.interp(sliceMaxMpc*100./LIGHT_SPEED, zVsDX, zVsDY)

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
#       ss0.000175 ~ SDSS DR9 mid 
baseResolution = float(numPart)/lbox/lbox/lbox # particles/Mpc^3
prevSubSample = -1
for thisSubSample in sorted(subSamples, reverse=True):

  keepFraction = float(thisSubSample) / baseResolution
  maxKeep = keepFraction * numPart
  minRadius = int(np.ceil(lbox/maxKeep**(1./3)))

  if args.script or args.all:
    print " Doing subsample", thisSubSample, "scripts"
    sys.stdout.flush()
    setName = prefix+"ss"+str(thisSubSample)
    if dataFormat == "multidark":
      subSampleToUse = 1.0
      fileToUse = prefix+"ss"+str(thisSubSample)+"_z"
      suffix = ".dat"
    elif dataFormat == "gadget":
      subSampleToUse = keepFraction
      #subSampleToUse = thisSubSample
      fileToUse = particleFileBase
      suffix = ""
    elif dataFormat == "lanl":
      subSampleToUse = keepFraction
      #subSampleToUse = thisSubSample
      fileToUse = particleFileBase
      suffix = ""
    elif dataFormat == "random":
      subSampleToUse = 1.0
      fileToUse = "ran.ss"+str(thisSubSample)+"_z"
      suffix = ".dat"

    writeScript(setName, fileToUse, dataFormat, 
                  scriptDir, catalogDir, fileNums, redshifts, 
                  numSubvolumes, numSlices, False, lbox, minRadius, omegaM, 
                  subsample=subSampleToUse, suffix=suffix)
    writeScript(setName, fileToUse, dataFormat, 
                  scriptDir, catalogDir, fileNums, redshifts, 
                  numSubvolumes, numSlices, True, lbox, minRadius, omegaM, 
                  subsample=subSampleToUse, suffix=suffix)
  
  if args.subsample or args.all:
    print " Doing subsample", thisSubSample
    sys.stdout.flush()

    for (iRedshift, redshift) in enumerate(redshifts):
      print "   redshift", redshift
      sys.stdout.flush()
  
      if dataFormat == "multidark":
        # reuse previous subamples in order to:
        #   - preserve unique IDs across multiple subsamples
        #   - reuse smaller files for faster processing
        if prevSubSample == -1:
          dataFile = catalogDir+"/"+particleFileBase+fileNums[iRedshift]
        else:
          sampleName = prefix+"ss"+str(prevSubSample)+"_z"+redshift
          dataFile = catalogDir+"/"+sampleName+".dat"
          keepFraction = float(thisSubSample) / float(prevSubSample)

        inFile = open(dataFile, 'r')

        sampleName = prefix+"ss"+str(thisSubSample)+"_z"+redshift
        outFile = open(catalogDir+"/"+sampleName+".dat", 'w')

        outFile.write("%f\n" %(lbox))
        outFile.write("%s\n" %(omegaM))
        outFile.write("%s\n" %(hubble))
        outFile.write("%s\n" %(redshift))
        outFile.write("%d\n" %(maxKeep))

        numKept = 0
        for (i,line) in enumerate(inFile):
          if (prevSubSample != -1 and i < 5): continue # skip header

          if np.random.uniform() > keepFraction: continue
          numKept += 1
          #if numKept > maxKeep: break

          if (prevSubSample == -1):
            line = line.split(',')
            x  = float(line[0])
            y  = float(line[1])
            z  = float(line[2])
            vz = float(line[3])
            vy = float(line[4])
            vx = float(line[5])
            uniqueID = i
            outFile.write("%d %e %e %e %e %e %e\n" %(uniqueID,x,y,z,vz,vy,vx))
          else:
            outFile.write(line)

        outFile.write("-99 -99 -99 -99 -99 -99 -99\n")
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

          outFile.write("%d %e %e %e 0. 0. 0.\n" % (i, x,y,z))

        outFile.write("-99 -99 -99 -99 -99 -99 -99\n")
        outFile.close()

  prevSubSample = thisSubSample

# -----------------------------------------------------------------------------
# now halos
if (args.script or args.all) and haloFileBase != "":

  for minHaloMass in minHaloMasses:
    print " Doing halo script", minHaloMass
    sys.stdout.flush()

    # estimate number of halos to get density
    if haloFileDummy == '':
      dataFile = catalogDir+haloFileBase+fileNums[0]
    else:
      dataFile = catalogDir+haloFileBase.replace(haloFileDummy, fileNums[0])

    inFile = open(dataFile, 'r')
    numPart = 0
    for (iLine, line) in enumerate(inFile): 
      if iLine < haloFileNumComLines: continue
      line = line.split(haloFileColSep)
      if minHaloMass == "none" or float(line[haloFileMCol]) > minHaloMass:
        numPart += 1
    inFile.close()

    minRadius = 2*int(np.ceil(lbox/numPart**(1./3.)))
  
    setName = prefix+"halos_min"+str(minHaloMass)
    writeScript(setName, prefix+"halos_min"+str(minHaloMass)+"_z", "multidark",
                scriptDir, catalogDir, fileNums, 
                redshifts, 
                numSubvolumes, numSlices, False, lbox, minRadius, omegaM)
    writeScript(setName, prefix+"halos_min"+str(minHaloMass)+"_z", "multidark",
                scriptDir, catalogDir, fileNums, 
                redshifts, 
                numSubvolumes, numSlices, True, lbox, minRadius, omegaM)

if (args.halos or args.all) and haloFileBase != "":
  print " Doing halos"
  sys.stdout.flush()

  for minHaloMass in minHaloMasses:
    print "  min halo mass = ", minHaloMass
    sys.stdout.flush()

    for (iRedshift, redshift) in enumerate(redshifts):
      print "   z = ", redshift
      sys.stdout.flush()

      if haloFileDummy == '':
        dataFile = catalogDir+haloFileBase+fileNums[iRedshift]
      else:
        dataFile = catalogDir+haloFileBase.replace(haloFileDummy, 
                                                   fileNums[iRedshift])
      inFile = open(dataFile, 'r')
      numPart = 0
      for (iLine, line) in enumerate(inFile): 
        if iLine < haloFileNumComLines: continue
        line = line.split(haloFileColSep)
        if minHaloMass == "none" or float(line[haloFileMCol]) > minHaloMass:
          numPart += 1
      inFile.close()

      sampleName = prefix+"halos_min"+str(minHaloMass)+"_z"+fileNums[iRedshift]
      outFile = open(catalogDir+"/"+sampleName+".dat", 'w')

      outFile.write("%f\n" %(lbox))
      outFile.write("%s\n" %(omegaM))
      outFile.write("%s\n" %(hubble))
      outFile.write("%s\n" %(redshift))
      outFile.write("%d\n" %(numPart))

      inFile = open(dataFile, 'r')
      for (iHalo,line) in enumerate(inFile):
        if iHalo < haloFileNumComLines: continue
        line = line.split(haloFileColSep)
        if minHaloMass == "none" or float(line[haloFileMCol]) > minHaloMass:
          x  = float(line[haloFileXCol])
          y  = float(line[haloFileYCol])
          z  = float(line[haloFileZCol])
          vz = float(line[haloFileVZCol])
          vy = float(line[haloFileVYCol])
          vx = float(line[haloFileVXCol])

          # write to output file
          outFile.write("%d %e %e %e %e %e %e\n" %(iHalo,x,y,z,vz,vy,vx))

      outFile.write("-99 -99 -99 -99 -99 -99 -99\n")
      inFile.close()
      outFile.close()

# -----------------------------------------------------------------------------
# now the HOD
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
GALAXY_DENSITY {galden} % computed automatically if M_min set, use for sanity
M1 {M1}   
sigma_logM {sigma_logM}
alpha {alpha}
M_cut {Mcut}

% simulation info
real_space_xi 1
HOD 1
populate_sim 1
HaloFile {haloFile}
HaloFileFormat {haloFileFormat}
RESOLUTION {numPartPerSide}
BOX_SIZE   {boxSize}

% output
root_filename {workDir}/hod
               """

if (args.script or args.all) and haloFileBase != "":
  print " Doing HOD scripts"
  sys.stdout.flush()
  for thisHod in hodParmList:
    print "   ", thisHod['name']
    setName = prefix+"hod_"+thisHod['name']
    writeScript(setName, prefix+"hod_"+thisHod['name']+"_z", "multidark",
                scriptDir, catalogDir, fileNums, redshifts, 
                 numSubvolumes, numSlices, False, lbox, 15, omegaM)
    writeScript(setName, prefix+"hod_"+thisHod['name']+"_z", "multidark",
                scriptDir, catalogDir, fileNums, redshifts, 
                 numSubvolumes, numSlices, True, lbox, 15, omegaM)

if (args.hod or args.all) and haloFileBase != "":
  print " Doing HOD"
  sys.stdout.flush()
  for thisHod in hodParmList:
    print "   ", thisHod['name']
    sys.stdout.flush()
    for (iRedshift, redshift) in enumerate(redshifts):
      print "  z = ", redshift
      sys.stdout.flush()

      parFileName = "./hod.par"
      parFile = open(parFileName, 'w')
      if haloFileDummy == '':
        haloFile = catalogDir+haloFileBase+fileNums[iRedshift]
      else:
        haloFile = catalogDir+haloFileBase.replace(haloFileDummy,
                                                 fileNums[iRedshift])
      parFile.write(parFileText.format(omegaM=omegaM,
                                     hubble=hubble,
                                     redshift=redshift,
                                     Mmin=thisHod['Mmin'],
                                     #Mmin=1.23e13,
                                     M1=thisHod['M1'],
                                     sigma_logM=thisHod['sigma_logM'],
                                     alpha=thisHod['alpha'],
                                     Mcut=thisHod['Mcut'],
                                     galden=thisHod['galDens'],
                                     #galden=0.000225,
                                     haloFile=haloFile,
                                     haloFileFormat=dataFormat,
                                     numPartPerSide=numPart**(1/3.),
                                     boxSize=lbox,
                                     workDir=catalogDir))
      parFile.close()

      os.system(hodPath+" "+parFileName+">& /dev/null")

      sampleName = getSampleName(prefix+"hod_"+thisHod['name'], redshift, False)
      outFileName = catalogDir+"/"+sampleName+".dat"
      os.system("mv %s/hod.mock %s" % (catalogDir, outFileName))
      os.system("rm %s/hod.*" % catalogDir)

print " Done!"
