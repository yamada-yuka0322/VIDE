#!/usr/bin/env python
#+
#   VIDE -- Void IDentification and Examination -- ./python_tools/pipeline_source/prepareCatalogs.in.py
#   Copyright (C) 2010-2013 Guilhem Lavaux
#   Copyright (C) 2011-2013 P. M. Sutter
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

# prepares input catalogs based on multidark simulations
#   (borrows heavily from generateMock, but doesn't hold much in memory)
# also creates necessary analyzeVoids input files

import numpy as np
import os
import sys
import void_python_tools as vp
import argparse
import imp
import subprocess

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
parser.add_argument('--parm', dest='parm',
                   default="",
                   help='path to parameter file')
args = parser.parse_args()


defaultsFile = "@CMAKE_BINARY_DIR@/python_tools/pipeline_source/defaults.py"
parms = imp.load_source("name", defaultsFile)
globals().update(vars(parms))

filename = args.parm
print " Loading parameters from", filename
if not os.access(filename, os.F_OK):
  print "  Cannot find parameter file %s!" % filename
  exit(-1)
parms = imp.load_source("name", filename)
globals().update(vars(parms))

#------------------------------------------------------------------------------
def getSDFTags(dataFile):
  idTag = "id"
  for line in open(dataFile):
    if "int64_t id" in line: 
      idTag = "id"
      break
    if "int64_t ident" in line: 
      idTag = "ident"
      break

  massTag = "mass"
  for line in open(dataFile):
    if "float m200b" in line: 
      massTag = "m200b"
      break
    if "float mass" in line: 
      massTag = "mass"
      break

  pidTag = ""
  iLine = 0
  for line in open(dataFile):
    iLine += 1
    if iLine > 100: break
    if "parent_id" in line: 
      pidTag = "parent_id"
      break
    if "pid" in line: 
      pidTag = "pid"
      break

  return idTag, massTag, pidTag

#------------------------------------------------------------------------------
def getSampleName(setName, redshift, useVel, iSlice=-1, iVol=-1):

  sampleName = setName

  sampleName += "_z" + redshift

  if iVol != -1: sampleName += "_d" + iVol

  return sampleName

#------------------------------------------------------------------------------
def getNickName(setName, sampleName):

  splitName = sampleName.split('_')

  nickName = datasetName  

  if "ss" in splitName[1]:
    nickName += " SS " + splitName[1].replace("ss","")
    #nickName = "Subsample = " + splitName[1].replace("ss","")
    if "pv" in splitName[2]: 
      nickName += ", z = " + splitName[3].replace("z","") + " (w/ PV)"
    else:
      nickName += ", z = " + splitName[2].replace("z","")
  elif "hod" in splitName[1]:
    nickName += " HOD " + splitName[2]
    if "pv" in splitName[3]: 
      nickName += ", z = " + splitName[4].replace("z","") + " (w/ PV)"
    else:
      nickName += ", z = " + splitName[3].replace("z","")
  elif "halos" in splitName[1]:
    if "none" in splitName[2]:
      nickName += " All Halos"
    else:
      nickName += " Halos > " + splitName[2].replace("min","").replace("+","")
    if "pv" in splitName[3]: 
      nickName += ", z = " + splitName[4].replace("z","") + " (w/ PV)"
    else:
      nickName += ", z = " + splitName[3].replace("z","")
  else:
    nickName = sampleName

  if "ran" in setName: nickName = "Random" + nickName

  return nickName


#------------------------------------------------------------------------------
# for given dataset parameters, outputs a script for use with analyzeVoids
def writeScript(setName, dataFileNameBase, dataFormat,
                scriptDir, catalogDir, fileNums, redshifts, numSubvolumes,
                numSlices, useVel, lbox, minRadius, omegaM, subsample=1.0, 
                suffix=".dat", dataFileNameList=None):


  if useVel: setName += "_pv"

  scriptFileName = scriptDir + "/" + setName + ".py"
  scriptFile = open(scriptFileName, 'w')

  header = """#!/usr/bin/env/python
import os
from void_python_tools.backend.classes import *

continueRun = {continueRun} # set to True to enable restarting aborted jobs
startCatalogStage = {startCatalogStage}
endCatalogStage   = {endCatalogStage}
               
startAPStage = {startAPStage}
endAPStage = {endAPStage}

regenerateFlag = False
#ZOBOV_PATH = "@CMAKE_BINARY_DIR@/c_tools/zobov2/"
ZOBOV_PATH = "@CMAKE_BINARY_DIR@/zobov/"
CTOOLS_PATH = "@CMAKE_BINARY_DIR@/c_tools/"
freshStack = True
errorBars = "CALCULATED"
numIncoherentRuns = 100
ranSeed = 101010
useComoving = False
bias = 1.16

dataPortions = {dataPortions}
dataSampleList = []
           """

  scriptFile.write(header.format(startCatalogStage=startCatalogStage,
                                 endCatalogStage=endCatalogStage,
                                 startAPStage=startAPStage,
                                 endAPStage=endAPStage,
                                 continueRun=continueRun,
                                 dataPortions=dataPortions))

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
                   nickName = "{nickName}",
                   dataType = "simulation",
                   zBoundary = ({zMin}, {zMax}),
                   zRange    = ({zMin}, {zMax}),
                   zBoundaryMpc = ({zMinMpc}, {zMaxMpc}),
                   omegaM    = {omegaM},
                   minVoidRadius = {minRadius},
                   profileBinSize = "auto",
                   includeInHubble = True,
                   partOfCombo = False,
                   {autoStack}
                   isCombo = False,
                   boxLen = {boxLen},
                   usePecVel = {usePecVel},
                   numSubvolumes = {numSubvolumes},
                   mySubvolume = "{mySubvolume}",
                   useLightCone = {useLightCone},
                   subsample = "{subsample}")
dataSampleList.append(newSample)
  """

  if stackMode == "fixed":
    stackInfo = """
# {zMin}, {zMax}, {minRadius}
newSample.addStack(0.0, 5.0, 5 , 10, False, False, rescaleMode="rv")
newSample.addStack(0.0, 5.0, 10, 15, False, False, rescaleMode="rv")
newSample.addStack(0.0, 5.0, 15, 20, False, False, rescaleMode="rv")
newSample.addStack(0.0, 5.0, 20, 25, False, False, rescaleMode="rv")
newSample.addStack(0.0, 5.0, 25, 30, False, False, rescaleMode="rv")
newSample.addStack(0.0, 5.0, 30, 35, False, False, rescaleMode="rv")
newSample.addStack(0.0, 5.0, 35, 40, False, False, rescaleMode="rv")
newSample.addStack(0.0, 5.0, 40, 45, False, False, rescaleMode="rv")
newSample.addStack(0.0, 5.0, 45, 50, False, False, rescaleMode="rv")
newSample.addStack(0.0, 5.0, 50, 55, False, False, rescaleMode="rv")
newSample.addStack(0.0, 5.0, 55, 60, False, False, rescaleMode="rv")
newSample.addStack(0.0, 5.0, 60, 65, False, False, rescaleMode="rv")
newSample.addStack(0.0, 5.0, 65, 70, False, False, rescaleMode="rv")
newSample.addStack(0.0, 5.0, 70, 75, False, False, rescaleMode="rv")
newSample.addStack(0.0, 5.0, 75, 80, False, False, rescaleMode="rv")
newSample.addStack(0.0, 5.0, 80, 85, False, False, rescaleMode="rv")
newSample.addStack(0.0, 5.0, 85, 90, False, False, rescaleMode="rv")
newSample.addStack(0.0, 5.0, 90, 95, False, False, rescaleMode="rv")
newSample.addStack(0.0, 5.0, 95, 100, False, False, rescaleMode="rv")
  """
  elif stackMode == "sample_thick":
    stackInfo = """
# {zMin}, {zMax}, {minRadius}
newSample.addStack({zMin}, {zMax}, 10, 20, True, False, rescaleMode="rmax")
newSample.addStack({zMin}, {zMax}, 20, 30, True, False, rescaleMode="rmax")
newSample.addStack({zMin}, {zMax}, 30, 40, True, False, rescaleMode="rmax")
newSample.addStack({zMin}, {zMax}, 40, 50, True, False, rescaleMode="rmax")
newSample.addStack({zMin}, {zMax}, 50, 60, True, False, rescaleMode="rmax")
newSample.addStack({zMin}, {zMax}, 60, 70, True, False, rescaleMode="rmax")
newSample.addStack({zMin}, {zMax}, 70, 80, True, False, rescaleMode="rmax")
newSample.addStack({zMin}, {zMax}, 80, 90, True, False, rescaleMode="rmax")
newSample.addStack({zMin}, {zMax}, 90, 100, True, False, rescaleMode="rmax")
  """

  elif stackMode == "log":
    stackInfo = ""
    rMin = 10
    rMax = 100
    rStart = rMin
    rEnd = rMin
    dlogR = 0.25
    while rEnd < rMax:
      rEnd = (1+0.5*dlogR)*rStart/(1-0.5*dlogR)
    
      stackInfo += """newSample.addStack({zMin}, {zMax}"""+ ", %g, %g, True, False, rescaleMode='rv')" % (rStart, rEnd)
      
      rStart = rEnd

  elif stackMode == "auto":
    stackInfo = """
newSample.addStack({zMin}, {zMax}, 2*{minRadius}  , 2*{minRadius}+2, True, False, rescaleMode="rv")
newSample.addStack({zMin}, {zMax}, 2*{minRadius}  , 2*{minRadius}+4, True, False, rescaleMode="rv")
newSample.addStack({zMin}, {zMax}, 2*{minRadius}+2, 2*{minRadius}+6, True, False, rescaleMode="rv")
newSample.addStack({zMin}, {zMax}, 2*{minRadius}+6, 2*{minRadius}+10, True, False, rescaleMode="rv")
newSample.addStack({zMin}, {zMax}, 2*{minRadius}+10, 2*{minRadius}+18, True, False, rescaleMode="rv")
newSample.addStack({zMin}, {zMax}, 2*{minRadius}+18, 2*{minRadius}+24, True, False, rescaleMode="rv")
               """
  else:
    stackInfo = """
# {zMin}, {zMax}, {minRadius}
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
      sliceMinMpc = "%0.2f" % sliceMinMpc
      sliceMaxMpc = "%0.2f" % sliceMaxMpc

      if (dataFileNameList != None):
        dataFileName = dataFileNameList[iFile]
      else:
        dataFileName = dataFileNameBase + fileNum + suffix

      for iX in xrange(numSubvolumes):
        for iY in xrange(numSubvolumes):

          mySubvolume = "%d%d" % (iX, iY)

          sampleName = getSampleName(setName, sliceMin, useVel,
                                     iSlice=iSlice, iVol=mySubvolume)
          nickName = getNickName(setName, sampleName)


          autoStack = ""
          if stackMode == "auto": autoStack = "autoNumInStack = 400,"
          scriptFile.write(sampleInfo.format(dataFile=dataFileName,
                                         dataFormat=dataFormat,
                                         dataUnit=dataUnit,
                                         sampleName=sampleName,
                                         nickName=nickName,
                                         zMin=sliceMin,
                                         zMax=sliceMax,
                                         zMinMpc=sliceMinMpc,
                                         zMaxMpc=sliceMaxMpc,
                                         omegaM=Om,
                                         boxLen=lbox,
                                         autoStack=autoStack,
                                         usePecVel=useVel,
                                         minRadius=minRadius,
                                         numSubvolumes=numSubvolumes,
                                         mySubvolume=mySubvolume,
                                         useLightCone=useLightCone,
                                         subsample=str(subsample).strip('[]')))

          for iAPSlice in xrange(numAPSlices):
            sliceWidth = float(sliceMax) - float(sliceMin)
            sliceAPMin = float(sliceMin) + iAPSlice*sliceWidth/numAPSlices
            sliceAPMax = float(sliceMin) + (iAPSlice+1)*sliceWidth/numAPSlices
            sliceAPMin = "%0.2f" % sliceAPMin
            sliceAPMax = "%0.2f" % sliceAPMax
            scriptFile.write(stackInfo.format(zMin=sliceAPMin,
                                              zMax=sliceAPMax,
                                              minRadius=minRadius))


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
subSamples = sorted(subSamples, reverse=True)
for iSubSample in xrange(len(subSamples)):

  subSampleList = subSamples[0:iSubSample+1]

  keepFractionList = []
  for subSample in subSampleList:
    if subSampleMode == 'absolute':
      keepFractionList.append(float(subSample) / baseResolution)
    elif subSampleMode == 'relative':
      keepFractionList.append(float(subSample))
    else:
      print "Unrecognized subSampleMode = ", subSampleMode
      exit(-1)
  thisSubSample = subSamples[iSubSample]
  maxKeep = keepFractionList[-1] * numPart
  minRadius = int(np.ceil(lbox/maxKeep**(1./3)))

  partFileList = []
  for (iRedshift, redshift) in enumerate(redshifts):
    if particleFileDummy == '':
      partFileList.append(particleFileBase+fileNums[iRedshift])
    else:
      partFileList.append(particleFileBase.replace(particleFileDummy, 
                                                   fileNums[iRedshift]))

  if args.script or args.all:
    print " Doing subsample", thisSubSample, "scripts"
    sys.stdout.flush()
    setName = prefix+"ss"+str(thisSubSample)

    if dataFormat == "random":
      subSampleToUse = 1.0
      fileToUse = "ran.ss"+str(thisSubSample)+"_z"
      suffix = ".dat"

      writeScript(setName, fileToUse, dataFormat, 
                  scriptDir, catalogDir, fileNums, redshifts, 
                  numSubvolumes, numSlices, False, lbox, minRadius, omegaM, 
                  subsample=subSampleToUse, suffix=suffix)
      if doPecVel:
        writeScript(setName, fileToUse, dataFormat, 
                    scriptDir, catalogDir, fileNums, redshifts, 
                    numSubvolumes, numSlices, True, lbox, minRadius, omegaM, 
                    subsample=subSampleToUse, suffix=suffix)
    else:
      if doSubSampling:
        # prepare scripts using our own format
        dataFormatToUse = "multidark"
        subSampleToUse = 1.0
        fileToUse = prefix+"ss"+str(thisSubSample)+"_z"
        partFileList = []
        suffix = ".dat"
        for (iRedshift, redshift) in enumerate(redshifts):
          sampleName = getSampleName(setName, redshift, False)
          outFileName = sampleName+".dat"
          partFileList.append(outFileName) 
      else:
        dataFormatToUse = dataFormat
        subSampleToUse = keepFractionList
        suffix = ""
        fileToUse = partFileList[0]
      writeScript(setName, fileToUse, dataFormatToUse, 
                  scriptDir, catalogDir, fileNums, redshifts, 
                  numSubvolumes, numSlices, False, lbox, minRadius, omegaM, 
                  subsample=subSampleToUse, suffix=suffix,
                  dataFileNameList=partFileList)
      if doPecVel:
        writeScript(setName, fileToUse, dataFormatToUse, 
                    scriptDir, catalogDir, fileNums, redshifts, 
                    numSubvolumes, numSlices, True, lbox, minRadius, omegaM, 
                    subsample=subSampleToUse, suffix=suffix,
                    dataFileNameList=partFileList)
       

  if (args.subsample or args.all) and doSubSampling:
    print " Doing subsample", thisSubSample
    sys.stdout.flush()

    for (iRedshift, redshift) in enumerate(redshifts):
      print "   redshift", redshift
      sys.stdout.flush()
  
      if dataFormat == "multidark" or dataFormat == "sdf":
        # reuse previous subamples in order to:
        #   - preserve unique IDs across multiple subsamples
        #   - reuse smaller files for faster processing
        if prevSubSample == -1:
          if particleFileDummy == '':
            dataFile = catalogDir+"/"+particleFileBase+fileNums[iRedshift]
          else:
            dataFile = particleFileBase.replace(particleFileDummy,
                                                fileNums[iRedshift])
            dataFile = catalogDir+"/"+dataFile
          keepFraction = float(thisSubSample) / baseResolution
        else:
          sampleName = prefix+"ss"+str(prevSubSample)+"_z"+redshift
          dataFile = catalogDir+"/"+sampleName+".dat"
          keepFraction = float(thisSubSample) / float(prevSubSample)

        if prevSubSample == -1 and dataFormat == "sdf":
          convertedFile = dataFile + "_temp"
          SDFcvt_PATH = "@CMAKE_BINARY_DIR@/external/libsdf/apps/SDFcvt/SDFcvt.x86_64"
          scale = 1./(1.+float(redshift))
          rescale_position = hubble/1000./scale
          shift = lbox/2.
          rescale_velocity = 3.08567802e16/3.1558149984e16

          command = "%s -a 200000 %s x y z vz vy vx mass | awk '{print $1*%g+%g, $2*%g+%g, $3*%g+%g, $4*%g, $5*%g, $6*%g, $7}' > %s" % (SDFcvt_PATH, dataFile,
                                     rescale_position,
                                     shift,
                                     rescale_position,
                                     shift,
                                     rescale_position,
                                     shift,
                                     rescale_velocity,
                                     rescale_velocity,
                                     rescale_velocity,
                                     convertedFile )
          #os.system(command)
          subprocess.call(command, shell=True)
          dataFile = convertedFile

        inFile = open(dataFile, 'r')

        sampleName = prefix+"ss"+str(thisSubSample)+"_z"+redshift
        outFile = open(catalogDir+"/"+sampleName+".dat", 'w')

        outFile.write("%f\n" %(lbox))
        outFile.write("%s\n" %(omegaM))
        outFile.write("%s\n" %(hubble))
        outFile.write("%s\n" %(redshift))
        outFile.write("%d\n" %(maxKeep))

        if dataFormat == "sdf":
          splitter = ' '
        if dataFormat == "multidark":
          splitter = ','

        numKept = 0
        for (i,line) in enumerate(inFile):
          if (prevSubSample != -1 and i < 5): continue # skip header

          if np.random.uniform() > keepFraction: continue
          numKept += 1
          #if numKept > maxKeep: break

          if (prevSubSample == -1):
            line = line.split(splitter)
            x  = float(line[0])
            y  = float(line[1])
            z  = float(line[2])
            vz = float(line[3])
            vy = float(line[4])
            vx = float(line[5])
            mass = float(line[6])
            uniqueID = i
            outFile.write("%d %e %e %e %e %e %e %e\n" %(uniqueID,x,y,z,
                                                     vz,vy,vx,mass))
          else:
            outFile.write(line)

        outFile.write("-99 -99 -99 -99 -99 -99 -99 -99\n")
        inFile.close()
        outFile.close()

        if prevSubSample == -1 and dataFormat == "sdf":
          os.unlink(dataFile)     
 
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

          outFile.write("%d %e %e %e 0. 0. 0. 0.\n" % (i, x,y,z))

        outFile.write("-99 -99 -99 -99 -99 -99 -99 -99\n")
        outFile.close()

  prevSubSample = thisSubSample

# -----------------------------------------------------------------------------
# now halos - filter on mass
if (args.script or args.all) and haloFileBase != "":

  for minHaloMass in minHaloMasses:
    print " Doing halo script", minHaloMass
    sys.stdout.flush()

    # estimate number of halos to get density
    # take from *highest* redshift to get lowest density
    if haloFileDummy == '':
      dataFile = catalogDir+haloFileBase+fileNums[-1]
    else:
      dataFile = catalogDir+haloFileBase.replace(haloFileDummy, 
                                                 fileNums[-1])
    numPart = 0
    if dataFormat == "sdf":
      SDFcvt_PATH = "@CMAKE_BINARY_DIR@/external/libsdf/apps/SDFcvt/SDFcvt.x86_64"
      idTag, massTag, pidTag = getSDFTags(dataFile)

      if minHaloMass == "none":
        command = "%s -a 200000 %s %s %s | awk '{if ($2==-1) print $1}' | wc" % (SDFcvt_PATH, dataFile, massTag, pidTag)
      else:
        command = "%s -a 200000 %s %s %s | awk '{if ($1>%g && $2==-1) print $1}' | wc" % (SDFcvt_PATH, dataFile, massTag, pidTag, minHaloMass)
      numPart = subprocess.check_output(command, shell=True)
      numPart = int(numPart.split()[0])
    else:
      inFile = open(dataFile, 'r')
      for (iHalo,line) in enumerate(inFile):
        if iHalo < haloFileNumComLines: continue
        line = line.split(haloFileColSep)
        if minHaloMass == "none" or float(line[haloFileMCol]) > minHaloMass:
          numPart += 1
      inFile.close()

    meanDen = lbox/numPart**(1./3.)
    if meanDen < 1: 
      minRadius = 4
    else:
      minRadius = int(np.ceil(meanDen))
  
    if minHaloMass != "none": 
      strMinHaloMass = "%.2e" % minHaloMass
    else:
      strMinHaloMass = "none"
 
    setName = prefix+"halos_min"+strMinHaloMass
    fileList = []
    for (iRedshift, redshift) in enumerate(redshifts):
      sampleName = getSampleName(setName, redshift, False)
      outFileName = sampleName+".dat"
      fileList.append(outFileName) 

    writeScript(setName, prefix+"halos_min"+strMinHaloMass+"_z", "multidark",
                scriptDir, catalogDir, fileNums, 
                redshifts, 
                numSubvolumes, numSlices, False, lbox, minRadius, omegaM,
                dataFileNameList = fileList)
   
    if doPecVel:
      writeScript(setName, prefix+"halos_min"+strMinHaloMass+"_z", 
                  "multidark",
                  scriptDir, catalogDir, fileNums, 
                  redshifts, 
                  numSubvolumes, numSlices, True, lbox, minRadius, omegaM,
                  dataFileNameList = fileList)

if (args.halos or args.all) and haloFileBase != "" and len(minHaloMasses) > 0:
  print " Doing halos - mass"
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
      if dataFormat == "sdf":
        for line in inFile:
          if "nhalos" in line:
            numPart = int(line.split()[3].strip(';'))
            break
          if "npart" in line:
            numPart = int(line.split()[3].strip(';'))
            break
        inFile.close()
      else:
        for (iLine, line) in enumerate(inFile):
          if iLine < haloFileNumComLines: continue
          line = line.split(haloFileColSep)
          if minHaloMass == "none" or float(line[haloFileMCol]) > minHaloMass:
            numPart += 1
        inFile.close()

      if minHaloMass != "none": 
        strMinHaloMass = "%.2e" % minHaloMass
      else:
        strMinHaloMass = "none"

      sampleName = prefix+"halos_min"+strMinHaloMass+"_z"+redshifts[iRedshift]
      outFileName = catalogDir+"/"+sampleName+".dat"
      outFile = open(outFileName, 'w')
      outFile.write("%f\n" %(lbox))
      outFile.write("%s\n" %(omegaM))
      outFile.write("%s\n" %(hubble))
      outFile.write("%s\n" %(redshift))
      outFile.write("%d\n" %(numPart))
      outFile.close()

      if dataFormat == "sdf":
        idTag, massTag, pidTag = getSDFTags(dataFile)

        SDFcvt_PATH = "@CMAKE_BINARY_DIR@/external/libsdf/apps/SDFcvt/SDFcvt.x86_64"
        if minHaloMass == "none":
          command = "%s -a 200000 %s %s %s x y z vz vy vx %s | awk '{if ($9==-1) print $2, $3, $4, $5, $6, $7, $8, $1}'>>%s" % (SDFcvt_PATH, dataFile, massTag, idTag, pidTag, outFileName )
        else:
          command = "%s -a 200000 %s %s %s x y z vz vy vx %s | awk '{if ($1>%g && $9==-1) print $2, $3, $4, $5, $6, $7, $8, $1}'>>%s" % (SDFcvt_PATH, dataFile, massTag, idTag, pidTag, minHaloMass, outFileName )
        #os.system(command)
        subprocess.call(command, shell=True)
        outFile = open(outFileName, 'a')
        outFile.write("-99 -99 -99 -99 -99 -99 -99 -99\n")
        outFile.close()
      else:
        outFile = open(outFileName, 'a')
        inFile = open(dataFile, 'r')
        for (iHalo,line) in enumerate(inFile):
          if iHalo < haloFileNumComLines: continue
          line = line.split(haloFileColSep)
          if minHaloMass == "none" or float(line[haloFileMCol]) > minHaloMass:
            x  = float(line[haloFileXCol]) * haloFilePosRescale
            y  = float(line[haloFileYCol]) * haloFilePosRescale
            z  = float(line[haloFileZCol]) * haloFilePosRescale
            vz = float(line[haloFileVZCol])
            vy = float(line[haloFileVYCol])
            vx = float(line[haloFileVXCol])
            mass = float(line[haloFileMCol])

            # write to output file
            outFile.write("%d %e %e %e %e %e %e %e\n" %(iHalo,x,y,z,
                                                        vz,vy,vx,mass))

        outFile.write("-99 -99 -99 -99 -99 -99 -99 -99\n")
        outFile.close()
        inFile.close()

# -----------------------------------------------------------------------------
# now halos - filter on density
if (args.script or args.all) and haloFileBase != "":

  for haloDen in haloDenList:
    print " Doing halo script", haloDen
    sys.stdout.flush()

    # estimate number of halos to get density
    if haloFileDummy == '':
      dataFile = catalogDir+haloFileBase+fileNums[0]
    else:
      dataFile = catalogDir+haloFileBase.replace(haloFileDummy, 
                                                 fileNums[0])

    numPart = haloDen * lbox**3
    meanDen = lbox/numPart**(1./3.)
    if meanDen < 1: 
      minRadius = 4
    else:
      minRadius = int(np.ceil(meanDen))
  
    setName = prefix+"halos_den"+str(haloDen)
    fileList = []
    for (iRedshift, redshift) in enumerate(redshifts):
      sampleName = getSampleName(setName, redshift, False)
      outFileName = sampleName+".dat"
      fileList.append(outFileName) 

    writeScript(setName, prefix+"halos_den"+str(haloDen)+"_z", "multidark",
                scriptDir, catalogDir, fileNums, 
                redshifts, 
                numSubvolumes, numSlices, False, lbox, minRadius, omegaM,
                dataFileNameList = fileList)
   
    if doPecVel:
      writeScript(setName, prefix+"halos_den"+str(haloDen)+"_z", 
                  "multidark",
                  scriptDir, catalogDir, fileNums, 
                  redshifts, 
                  numSubvolumes, numSlices, True, lbox, minRadius, omegaM,
                  dataFileNameList = fileList)

if (args.halos or args.all) and haloFileBase != "":
  print " Doing halos - density"
  sys.stdout.flush()

  for haloDen in haloDenList:
    print "  halo density = ", haloDen
    sys.stdout.flush()

    for (iRedshift, redshift) in enumerate(redshifts):
      print "   z = ", redshift
      sys.stdout.flush()

      if haloFileDummy == '':
        dataFile = catalogDir+haloFileBase+fileNums[iRedshift]
      else:
        dataFile = catalogDir+haloFileBase.replace(haloFileDummy, 
                                                   fileNums[iRedshift])
      # check the actual density of halos
      inFile = open(dataFile, 'r')
      numPart = 0
      if dataFormat == "sdf":
        for line in inFile:
          if "nhalos" in line:
            numPart = int(line.split()[3].strip(';'))
            break
          if "npart" in line:
            numPart = int(line.split()[3].strip(';'))
            break
        inFile.close()
      else:
        for (iLine, line) in enumerate(inFile):
          if iLine < haloFileNumComLines: continue
          numPart += 1
        inFile.close()

      numPartExpect = int(np.ceil(haloDen * lbox**3))
      if numPart < numPartExpect:
        print " ERROR: not enough halos to support that density! Maximum is %g" % (1.*numPart / lbox**3)
        exit(-1)


      sampleName = prefix+"halos_den"+str(haloDen)+"_z"+redshifts[iRedshift]
      outFileName = catalogDir+"/"+sampleName+".dat"
      outFile = open(outFileName, 'w')
      outFile.write("%f\n" %(lbox))
      outFile.write("%s\n" %(omegaM))
      outFile.write("%s\n" %(hubble))
      outFile.write("%s\n" %(redshift))
      outFile.write("%d\n" %(numPartExpect))
      outFile.close()

      if dataFormat == "sdf":
        idTag, massTag, pidTag = getSDFTags(dataFile)
      
        tempFile = catalogDir+"/temp_"+sampleName+".dat" 
        SDFcvt_PATH = "@CMAKE_BINARY_DIR@/external/libsdf/apps/SDFcvt/SDFcvt.x86_64"
        command = "%s -a 200000 %s %s %s x y z vz vy vx %s | awk '{if ($9==-1) print $2, $3, $4, $5, $6, $7, $8, $1}'>>%s" % (SDFcvt_PATH, dataFile, massTag, idTag, pidTag, tempFile )
        #os.system(command)
        subprocess.call(command, shell=True)

        # do the check again since now we've filtered out subhalos
        numPart = 0
        fileTemp = open(tempFile, 'r')
        for line in fileTemp:
          numPart += 1
        fileTemp.close()

        actualDen = 1.*numPart / lbox**3
        keepFraction = haloDen / actualDen
        if numPart < numPartExpect:
          print " ERROR: not enough galaxies to support that density! Maximum is %g" % (1.*numPart / lbox**3)
          exit(-1)

        numKept = 0
        inFile = open(tempFile, 'r')
        outFile = open(outFileName, 'a')
        for (iLine,line) in enumerate(inFile):
          if np.random.uniform() > keepFraction: continue
          outFile.write(line)
          numKept += 1
        inFile.close()

        outFile.write("-99 -99 -99 -99 -99 -99 -99 -99\n")
        outFile.close()
      else:
        actualDen = 1.*numPart / lbox**3
        keepFraction = haloDen / actualDen

        inFile = open(dataFile, 'r')
        outFile = open(outFileName, 'a')
        for (iHalo,line) in enumerate(inFile):
          if iHalo < haloFileNumComLines: continue
          if np.random.uniform() > keepFraction: continue
          line = line.split(haloFileColSep)
          x  = float(line[haloFileXCol]) * haloFilePosRescale
          y  = float(line[haloFileYCol]) * haloFilePosRescale
          z  = float(line[haloFileZCol]) * haloFilePosRescale
          vz = float(line[haloFileVZCol])
          vy = float(line[haloFileVYCol])
          vx = float(line[haloFileVXCol])
          mass = float(line[haloFileMCol])

          outFile.write("%d %e %e %e %e %e %e %e\n" %(iHalo,x,y,z,
                                                      vz,vy,vx,mass))
        outFile.write("-99 -99 -99 -99 -99 -99 -99 -99\n")
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
root_filename {workDir}/hod_{sampleName}
               """

if (args.script or args.all) and haloFileBase != "":
  print " Doing HOD scripts"
  sys.stdout.flush()

  for thisHod in hodParmList:
    fileList = []
    for (iRedshift, redshift) in enumerate(redshifts):
      sampleName = getSampleName(prefix+"hod_"+thisHod['name'], redshift, False)
      outFileName = sampleName+".dat"
      fileList.append(outFileName) 

    print "  ", thisHod['name']

    # estimate number of halos to get density
    numPart = thisHod['galDensFinal'] * lbox**3
    meanDen = lbox/numPart**(1./3.)
    if meanDen < 1: 
      minRadius = 4
    else:
      minRadius = int(np.ceil(meanDen))
  
    setName = prefix+"hod_"+thisHod['name']
    writeScript(setName, prefix+"hod_"+thisHod['name']+"_z", "multidark",
                scriptDir, catalogDir, fileNums, redshifts, 
                numSubvolumes, numSlices, False, lbox, minRadius, omegaM,
                dataFileNameList = fileList)
    if doPecVel:
      writeScript(setName, prefix+"hod_"+thisHod['name']+"_z", "multidark",
                  scriptDir, catalogDir, fileNums, redshifts, 
                  numSubvolumes, numSlices, True, lbox, minRadius, omegaM,
                  dataFileNameList = fileList)

if (args.hod or args.all) and haloFileBase != "":
  print " Doing HOD"
  sys.stdout.flush()
  for (iRedshift, redshift) in enumerate(redshifts):
    print "  z = ", redshift
    sys.stdout.flush()

    if haloFileDummy == '':
      haloFile = catalogDir+haloFileBase+fileNums[iRedshift]
    else:
      haloFile = catalogDir+haloFileBase.replace(haloFileDummy,
                                                 fileNums[iRedshift])

    if dataFormat == "sdf":
      inFile = haloFile
      outFile = haloFile+"_temp"

      idTag, massTag, pidTag = getSDFTags(inFile)

      SDFcvt_PATH = "@CMAKE_BINARY_DIR@/external/libsdf/apps/SDFcvt/SDFcvt.x86_64"
      command = "%s -a 200000 %s %s x y z vx vy vz %s | awk '{if ($8 ==-1) print $1, $2, $3, $4, $5, $6, $7}'>>%s" % (SDFcvt_PATH, inFile, massTag, pidTag, outFile)
      #os.system(command)
      subprocess.call(command, shell=True)
      haloFile = outFile

    for thisHod in hodParmList:
      print "   ", thisHod['name']
      sys.stdout.flush()

      sampleName = getSampleName(prefix+"hod_"+thisHod['name'], redshift, False)

      parFileName = "./hod.par"
      parFile = open(parFileName, 'w')
      parFile.write(parFileText.format(omegaM=omegaM,
                                     hubble=hubble,
                                     redshift=redshift,
                                     Mmin=thisHod['Mmin'],
                                     M1=thisHod['M1'],
                                     sigma_logM=thisHod['sigma_logM'],
                                     alpha=thisHod['alpha'],
                                     Mcut=thisHod['Mcut'],
                                     #galden=-1,
                                     galden=thisHod['galDens'],
                                     haloFile=haloFile,
                                     haloFileFormat=dataFormat,
                                     numPartPerSide=numPart**(1/3.),
                                     boxSize=lbox,
                                     workDir=catalogDir,
                                     sampleName=sampleName))
      parFile.close()

      tempFile = "./hod.out_"+sampleName
      output = open(tempFile, 'w')
      HOD_PATH = "@CMAKE_BINARY_DIR@/c_tools/hod/hod"
      
      #os.system(HOD_PATH+" "+parFileName+">& " + tempFile)
      subprocess.call(HOD_PATH+" "+parFileName, stdout=output, stderr=output, 
                      shell=True)
      output.close()
  
      hodWorked = False 
      for line in open(tempFile):
        if "MLO" in line:
          print "     (minimum halo mass = ", line.split()[1], ")"
          hodWorked = True 
          break

      if hodWorked:
        os.unlink(tempFile)
      else:
        print "HOD Failed! Log follows:"
        for line in open(tempFile):
          print line
        exit(-1)

      # now randomly subsample the galaxies to get desired density
      inFile = open(catalogDir+"/hod_"+sampleName+".mock")
      numPartExpect = thisHod['galDensFinal'] * lbox**3
      numPartActual = 0
      actualDen = 0
      for (iLine,line) in enumerate(inFile):
        numPartActual += 1
      inFile.close()

      if numPartActual < numPartExpect:
        print " ERROR: not enough galaxies to support that density! Maximum is %g" % (1.*numPartActual / lbox**3)
        exit(-1)

      actualDen = 1.*numPartActual / lbox**3
      keepFraction = float(thisHod['galDensFinal']) / actualDen

      inFile = open(catalogDir+"/hod_"+sampleName+".mock")
      outFile = open(catalogDir+"/"+sampleName+".dat", 'w')

      #outFile.write("%f\n" %(lbox))
      #outFile.write("%s\n" %(omegaM))
      #outFile.write("%s\n" %(hubble))
      #outFile.write("%s\n" %(redshift))
      #outFile.write("%d\n" %(numPartExpect))

      numKept = 0
      for (iLine,line) in enumerate(inFile):
        if iLine < 5: 
          outFile.write(line)
          continue
        if np.random.uniform() > keepFraction: continue
        outFile.write(line) 
        numKept += 1

      inFile.close()
      os.unlink(catalogDir+"/hod_"+sampleName+".mock") 

      outFile.write("-99 -99 -99 -99 -99 -99 -99 -99\n")
      outFile.close()
      

      #os.system("mv %s/hod_%s.mock %s" % (catalogDir, sampleName, outFileName))
      #os.system("rm %s/hod.*" % catalogDir)
      os.system("rm ./hod.par")
      os.system("rm ./hod-usedvalues")

    if dataFormat == "sdf": os.system("rm %s" % haloFile)

print " Done!"
