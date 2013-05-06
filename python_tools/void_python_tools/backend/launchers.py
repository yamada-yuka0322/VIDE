#+
#   VIDE -- Void IDEntification pipeline -- ./python_tools/void_python_tools/backend/launchers.py
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
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# routines which communicate with individual data analysis portions - 
#   they make the analyzeVoids.py script easier to read

import os
import glob
import classes
import numpy as np
import numpy.ma as ma
import os
import shutil
import glob
import subprocess
import sys
from   pylab import figure
from   netCDF4 import Dataset
from void_python_tools.backend.classes import *
import void_python_tools.apTools as vp

NetCDFFile = Dataset
ncFloat = 'f8' # Double precision

# -----------------------------------------------------------------------------
def launchGenerate(sample, binPath, workDir=None, inputDataDir=None, 
                   zobovDir=None, figDir=None, logFile=None, useLCDM=False,
                   continueRun=None,regenerate=False):

  if sample.dataType == "observation":
    sampleName = sample.fullName

    if regenerate:
      inputParameterFlag = "inputParameter " + zobovDir+"/zobov_slice_"+sampleName+".par"
      outputFile = zobovDir+"/regenerated_zobov_slice_" + sampleName
    else:
      inputParameterFlag = ""
      outputFile = zobovDir+"/zobov_slice_" + sampleName

    if sample.dataFile == "":
      datafile = inputDataDir+"/"+sampleName
    else:
      datafile = inputDataDir+"/"+sample.dataFile

    maskFile = sample.maskFile

    if useLCDM:
      useLCDMFlag = "useLCDM"
    else:
      useLCDMFlag = ""

    conf="""
      catalog %s
      mask %s
      output %s
      params %s
      zMin %g
      zMax %g
      density_fake %g
      %s
      %s
      """ % (datafile, maskFile, outputFile,
             zobovDir+"/zobov_slice_"+sampleName+".par",
             sample.zBoundary[0], sample.zBoundary[1], sample.fakeDensity,
             useLCDMFlag, inputParameterFlag)

    parmFile = os.getcwd()+"/generate_"+sample.fullName+".par"

    if regenerate or not (continueRun and jobSuccessful(logFile, "Done!\n")):
      file(parmFile, mode="w").write(conf)
      cmd = "%s --configFile=%s &> %s" % (binPath,parmFile,logFile)
      os.system(cmd)
      if jobSuccessful(logFile, "Done!\n"):
        print "done"
      else:
        print "FAILED!"
        exit(-1)

    else:
      print "already done!"

    if os.access(parmFile, os.F_OK): os.unlink(parmFile)

    if os.access("contour_map.fits", os.F_OK):
      os.system("mv %s %s" % ("contour_map.fits", zobovDir))
      os.system("mv %s %s" % ("mask_map.fits", zobovDir))

    if os.access("comoving_distance.txt", os.F_OK):
      os.system("mv %s %s" % ("comoving_distance.txt", zobovDir))

    if os.access("mask_index.txt", os.F_OK):
      os.system("mv %s %s" % ("mask_index.txt", zobovDir))
      os.system("mv %s %s" % ("total_particles.txt", zobovDir))
      #os.system("mv %s %s" % ("sample_info.txt", zobovDir))

    if os.access("galaxies.txt", os.F_OK):
      os.system("mv %s %s" % ("galaxies.txt", zobovDir))
      os.system("mv %s %s" % ("mock_galaxies.txt", zobovDir))
      os.system("mv %s %s" % ("mock_boundary.txt", zobovDir))
      os.system("mv %s %s" % ("mock_sphere.txt", zobovDir))

  else: # simulation
    sampleName = sample.fullName

    datafile = inputDataDir+"/"+sample.dataFile

    # check if the final subsampling is done
    lastSample = sample.subsample.split(', ')[-1]
    doneLine = "Done! %5.2e\n" % float(lastSample)
    if (continueRun and jobSuccessful(logFile, doneLine)): 
      print "already done!"
      return

    prevSubSample = -1
    firstSubSample = -1
    for thisSubSample in sample.subsample.split(', '):

      if prevSubSample == -1:
        inputParameterFlag = ""
        outputFile = zobovDir+"/zobov_slice_" + sampleName + "_ss" + \
                     thisSubSample
        keepFraction = float(thisSubSample)
        subSampleLine = "subsample %g" % keepFraction
        resubSampleLine = ""
        firstSubSample = keepFraction
      else:
        inputParameterFlag = "inputParameter " + zobovDir+"/zobov_slice_"+\
                             sampleName+"_ss"+prevSubSample+".par"
        outputFile = zobovDir+"/zobov_slice_" + sampleName + "_ss" + \
                     thisSubSample
        keepFraction = float(thisSubSample)/float(prevSubSample)
        subSampleLine = "subsample %s" % firstSubSample
        resubSampleLine = "resubsample %g" % keepFraction

      includePecVelString = ""
      if sample.usePecVel: includePecVelString = "peculiarVelocities"

      useLightConeString = ""
      if sample.useLightCone: useLightConeString = "cosmo"

      if sample.dataFormat == "multidark" or sample.dataFormat == "random":
        dataFileLine = "multidark " + datafile
      elif sample.dataFormat == "gadget":
        dataFileLine = "gadget " + datafile
      elif sample.dataFormat == "sdf":
        dataFileLine = "sdf " + datafile
      else:
        raise ValueError("unknown dataFormat '%s'" % sample.dataFormat)
    
      iX = float(sample.mySubvolume[0])
      iY = float(sample.mySubvolume[1])

      xMin = iX/sample.numSubvolumes * sample.boxLen
      yMin = iY/sample.numSubvolumes * sample.boxLen
      xMax = (iX+1)/sample.numSubvolumes * sample.boxLen
      yMax = (iY+1)/sample.numSubvolumes * sample.boxLen

      conf="""
      %s
      output %s
      outputParameter %s
      %s
      %s
      gadgetUnit %g
      rangeX_min %g
      rangeX_max %g
      rangeY_min %g
      rangeY_max %g
      rangeZ_min %g
      rangeZ_max %g
      %s
      %s
      %s
      joggleParticles
      """ % (dataFileLine, outputFile,
             outputFile+".par",
             includePecVelString,
             useLightConeString,
             sample.dataUnit,
             xMin, xMax, yMin, yMax,
             sample.zBoundaryMpc[0], sample.zBoundaryMpc[1],
             subSampleLine,resubSampleLine,inputParameterFlag)

      parmFile = os.getcwd()+"/generate_"+sample.fullName+".par"

      file(parmFile, mode="w").write(conf)

      if (prevSubSample == -1):
        cmd = "%s --configFile=%s &> %s" % (binPath,parmFile,logFile)
      else:
        cmd = "%s --configFile=%s &>> %s" % (binPath,parmFile,logFile)
      os.system(cmd)
 
      # remove intermediate files
      if (prevSubSample != -1):
       os.unlink(zobovDir+"/zobov_slice_"+sampleName+"_ss"+prevSubSample+".par")
       os.unlink(zobovDir+"/zobov_slice_"+sampleName+"_ss"+prevSubSample)

      doneLine = "Done! %5.2e\n" % keepFraction
      if not jobSuccessful(logFile, doneLine):
        print "FAILED!"
        exit(-1)

      prevSubSample = thisSubSample

    if jobSuccessful(logFile, doneLine): print "done"
 
    # place the final subsample     
    os.system("mv %s %s" % (zobovDir+"/zobov_slice_"+sampleName+"_ss"+\
              prevSubSample, zobovDir+"/zobov_slice_"+sampleName))
    os.system("mv %s %s" % (zobovDir+"/zobov_slice_"+sampleName+"_ss"+\
              prevSubSample+".par", zobovDir+"/zobov_slice_"+sampleName+".par"))

    if os.access("comoving_distance.txt", os.F_OK):
      os.system("mv %s %s" % ("comoving_distance.txt", zobovDir))

    if os.access(parmFile, os.F_OK):
      os.unlink(parmFile)

# -----------------------------------------------------------------------------
def launchZobov(sample, binPath, zobovDir=None, logDir=None, continueRun=None,
                 numZobovDivisions=None, numZobovThreads=None):

  sampleName = sample.fullName

  datafile = zobovDir+"zobov_slice_"+sampleName

  logFile = logDir+"/zobov_"+sampleName+".out"

  vozScript = "./scr_"+sampleName

  if os.access(vozScript, os.F_OK):
    os.unlink(vozScript)

  if sample.dataType == "observation":
    maskIndex = open(zobovDir+"/mask_index.txt", "r").read()
    totalPart = open(zobovDir+"/total_particles.txt", "r").read()
    maxDen = 0.2*float(maskIndex)/float(totalPart)
  else:
    maskIndex = -1
    maxDen = 0.2

  if not (continueRun and jobSuccessful(logFile, "Done!\n")):
    for fileName in glob.glob(zobovDir+"/part._"+sampleName+".*"):
      os.unlink(fileName)

    if os.access(zobovDir+"/adj_"+sampleName+".dat", os.F_OK):
      os.unlink(zobovDir+"adj_"+sampleName+".dat")

    if os.access(zobovDir+"/voidDesc_"+sampleName+".out", os.F_OK):
      os.unlink(zobovDir+"/voidDesc_"+sampleName+".out")

    cmd = "%s/vozinit %s 0.1 1.0 %g %s %g %s %s %s &> %s" % \
          (binPath, datafile,  numZobovDivisions, \
                      "_"+sampleName, numZobovThreads, \
                      binPath, zobovDir, maskIndex, logFile)
    os.system(cmd)

    cmd = "./%s >> %s 2>&1" % (vozScript, logFile)
    os.system(cmd)

    cmd = "%s/../c_tools/zobov2/jozov2 %s %s %s %s %s %g %s >> %s 2>&1" % \
          (binPath, \
           zobovDir+"/adj_"+sampleName+".dat", \
           zobovDir+"/vol_"+sampleName+".dat", \
           zobovDir+"/voidPart_"+sampleName+".dat", \
           zobovDir+"/voidZone_"+sampleName+".dat", \
           zobovDir+"/voidDesc_"+sampleName+".out", \
           maxDen, maskIndex, logFile)
    os.system(cmd)

    if jobSuccessful(logFile, "Done!\n"):
      print "done"
    else:
      print "FAILED!"
      exit(-1)

  else:
    print "already done!"

  if os.access(vozScript, os.F_OK):
    os.unlink(vozScript)
 
# -----------------------------------------------------------------------------
def launchPrune(sample, binPath, 
                summaryFile=None, logFile=None, zobovDir=None, 
                continueRun=None):

  sampleName = sample.fullName

  numVoids = sum(1 for line in \
                 open(zobovDir+"/voidDesc_"+sampleName+".out"))
  numVoids -= 2

  if sample.dataType == "observation":
    mockIndex = open(zobovDir+"/mask_index.txt", "r").read()
    totalPart = open(zobovDir+"/total_particles.txt", "r").read()
    maxDen = 0.2*float(mockIndex)/float(totalPart)
    observationLine = " --isObservation"
    #periodicLine = " --periodic=''"
  else:
    mockIndex = -1
    maxDen = 0.2
    observationLine = ""

    #periodicLine = " --periodic='"
    #if sample.numSubvolumes == 1: periodicLine += "xy"
    #if np.abs(sample.zBoundaryMpc[1]  - sample.zBoundaryMpc[0] - \
    #   sample.boxLen) <= 1.e-1: 
    #  periodicLine += "z"
    #periodicLine += "' "

  periodicLine = " --periodic='" + getPeriodic(sample) + "'"

  if not (continueRun and (jobSuccessful(logFile, "NetCDF: Not a valid ID\n") \
          or jobSuccessful(logFile, "Done!\n"))):
    cmd = binPath
    cmd += " --partFile=" + zobovDir+"/zobov_slice_"+str(sampleName)
    cmd += " --voidDesc=" + zobovDir+"/voidDesc_"+str(sampleName)+".out"
    cmd += " --void2Zone="+zobovDir+"/voidZone_"+str(sampleName)+".dat"
    cmd += " --zone2Part=" + zobovDir+"/voidPart_"+str(sampleName)+".dat"
    cmd += " --partVol=" + zobovDir+"/vol_"+str(sampleName)+".dat"
    cmd += " --extraInfo=" + zobovDir+"/zobov_slice_"+str(sampleName)+\
           ".par"
    cmd += " --tolerance=1.0"
    cmd += " --mockIndex=" + str(mockIndex)
    cmd += " --maxCentralDen=" + str(maxDen)
    cmd += " --zMin=" + str(sample.zRange[0])
    cmd += " --zMax=" + str(sample.zRange[1])
    cmd += " --rMin=" + str(sample.minVoidRadius)
    cmd += " --numVoids=" + str(numVoids)
    cmd += observationLine
    cmd += periodicLine
    cmd += " --outputDir=" + zobovDir
    cmd += " --sampleName=" + str(sampleName)
    cmd += " &> " + logFile
    f=file("run_prune.sh",mode="w")
    f.write(cmd)
    f.close()
    os.system(cmd)

    if jobSuccessful(logFile, "NetCDF: Not a valid ID\n") or \
       jobSuccessful(logFile, "Done!\n"):
      print "done"
    else:
      print "FAILED!"
      #exit(-1)

  else:
    print "already done!"

 
# -----------------------------------------------------------------------------
def launchVoidOverlap(sample1, sample2, sample1Dir, sample2Dir, 
                      binPath, thisDataPortion=None, 
                      logFile=None, workDir=None,
                      continueRun=None, outputFile=None, 
                      matchMethod=None):

  sampleName1 = sample1.fullName
  sampleName2 = sample2.fullName

  periodicLine = " --periodic='" + getPeriodic(sample1) + "' "
  #if sample1.dataType != "observation":
  #  if sample1.numSubvolumes == 1: periodicLine += "xy"
  #  if np.abs(sample1.zBoundaryMpc[1]  - sample1.zBoundaryMpc[0] - sample1.boxLen) <= 1.e-1: 
  #    periodicLine += "z"
  #periodicLine += "' "

  if not (continueRun and jobSuccessful(logFile, "Done!\n")):
    cmd = binPath
    cmd += " --partFile1=" + sample1Dir+"/zobov_slice_" + \
           str(sampleName1)
    cmd += " --volFile1=" + sample1Dir+"/vol_" + \
           str(sampleName1)+".dat"
    cmd += " --voidFile1=" + sample1Dir+"/voidDesc_" + \
           thisDataPortion+"_"+str(sampleName1)+".out"
    cmd += " --infoFile1=" + sample1Dir+"/zobov_slice_" + \
           str(sampleName1)+".par"
    cmd += " --centerFile1=" + sample1Dir + \
           "/barycenters_"+thisDataPortion+"_"+str(sampleName1)+".out"
    cmd += " --shapeFile1=" + sample1Dir + \
           "/shapes_"+thisDataPortion+"_"+str(sampleName1)+".out"
    cmd += " --zoneFile1=" + sample1Dir+"/voidZone_" + \
           str(sampleName1)+".dat"
    cmd += " --zonePartFile1=" + sample1Dir+"/voidPart_" + \
           str(sampleName1)+".dat"

    cmd += " --partFile2=" + sample2Dir+"/zobov_slice_" + \
           str(sampleName2)
    cmd += " --volFile2=" + sample2Dir+"/vol_" + \
           str(sampleName2)+".dat"
    cmd += " --voidFile2=" + sample2Dir+"/voidDesc_" + \
           thisDataPortion+"_"+str(sampleName2)+".out"
    cmd += " --infoFile2=" + sample2Dir+"/zobov_slice_" + \
           str(sampleName2)+".par"
    cmd += " --centerFile2=" + sample2Dir + \
           "/barycenters_"+thisDataPortion+"_"+str(sampleName2)+".out"
    cmd += " --shapeFile2=" + sample1Dir + \
           "/shapes_"+thisDataPortion+"_"+str(sampleName1)+".out"
    cmd += " --zoneFile2=" + sample2Dir+"/voidZone_" + \
           str(sampleName2)+".dat"
    cmd += " --zonePartFile2=" + sample2Dir+"/voidPart_" + \
           str(sampleName2)+".dat"

    if matchMethod == "useID": cmd += " --useID"
    cmd += periodicLine
    cmd += " --outfile=" + outputFile
    cmd += " &> " + logFile
    open("temp.par",'w').write(cmd)
    os.system(cmd)

    if jobSuccessful(logFile, "Done!\n"):
      print "done"
    else:
      print "FAILED!"
      exit(-1)

  else:
    print "already done!"


# -----------------------------------------------------------------------------
def launchStack(sample, stack, binPath, thisDataPortion=None, logDir=None,
                voidDir=None, freshStack=True, runSuffix=None,
                zobovDir=None,
                INCOHERENT=False, ranSeed=None, summaryFile=None, 
                continueRun=None, dataType=None, prefixRun="",
                idList=None, rescaleOverride=None):

  sampleName = sample.fullName

  if idList != None:
    customLine = "selected"
  else:
    customLine = ""

  runSuffix = getStackSuffix(stack.zMin, stack.zMax, stack.rMin,
                             stack.rMax, thisDataPortion, 
                             customLine=customLine)

  logFile = logDir+("/%sstack_"%prefixRun)+sampleName+"_"+runSuffix+".out"
 
  treeFile = voidDir+"/tree.data"

  if (freshStack) and os.access(treeFile, os.F_OK):
    os.unlink(treeFile)

  centralRadius = stack.rMin * 0.25

  # restrict to relavent ranges of sample
  zMin = max(sample.zRange[0],stack.zMin) * 3000
  zMax = min(sample.zRange[1],stack.zMax) * 3000

  if dataType == "observation":
    obsFlag = "observation"
  else:
    obsFlag = ""

  Rextracut = stack.rMin*3 + 1
  Rcircular = stack.rMin*3 + 2

  if dataType == "observation":
    maskIndex = open(zobovDir+"/mask_index.txt", "r").read()
    totalPart = open(zobovDir+"/total_particles.txt", "r").read()
    maxDen = 0.2*float(maskIndex)/float(totalPart)
  else:
    maskIndex = 999999999
    maxDen = 0.2

  if INCOHERENT:
    nullTestFlag = "INCOHERENT"
  else:
    nullTestFlag = ""

  if rescaleOverride == None: rescaleOverride = stack.rescaleMode
  if rescaleOverride == "rmax":
    rescaleFlag = "rescale"
  else:
    rescaleFlag = ""

  idListFile = "idlist.temp"
  if idList != None:
    idListFlag = "idList " + idListFile
    idFile = open(idListFile, 'w')
    for id in idList: idFile.write(str(id)+"\n")
    idFile.close()
    rMinToUse = 0.
    rMaxToUse = 100000.
    centralRadius = rMaxToUse
  else:
    idListFlag = ""
    rMinToUse = stack.rMin
    rMaxToUse = stack.rMax

  periodicLine = "periodic='" + getPeriodic(sample) + "'"
  #if sample.dataType == "observation":
  #  periodicLine = "periodic=''"
  #else:
  #  periodicLine = "periodic='"
  #  if sample.numSubvolumes == 1: periodicLine += "xy"
  #  if np.abs(sample.zBoundaryMpc[1]  - sample.zBoundaryMpc[0] - \
  #     sample.boxLen) <= 1.e-1:
  #    periodicLine += "z"
  #  periodicLine += "' "

  conf="""
  desc %s
  partzone %s
  zonevoid %s
  volumefile %s
  Rmin %g
  Rmax %g
  particles %s
  extraInfo %s
  densityThreshold %g
  centralRadius %g
  edgeAvoidance %g
  Circular %g
  Zmin %g
  Zmax %g
  %s
  %s
  ranSeed %d
  dataPortion %s
  barycenters %s
  boundaryDistances %s
  %s
  %s
  %s
  """ % \
  (zobovDir+"/voidDesc_"+thisDataPortion+"_"+sampleName+".out",
   zobovDir+"/voidPart_"+sampleName+".dat",
   zobovDir+"/voidZone_"+sampleName+".dat",
   zobovDir+"/vol_"+sampleName+".dat",
   rMinToUse,
   rMaxToUse,
   zobovDir+("/%szobov_slice_"%prefixRun)+sampleName,
   zobovDir+"/zobov_slice_"+sampleName+".par",
   maxDen,
   centralRadius,
   Rextracut,
   Rcircular,
   zMin,
   zMax,
   obsFlag,
   nullTestFlag,
   ranSeed,
   thisDataPortion,
   zobovDir+"/barycenters_"+thisDataPortion+"_"+sampleName+".out",
   zobovDir+"/boundaryDistances_"+thisDataPortion+"_"+sampleName+".out",
   rescaleFlag,
   idListFlag,
   periodicLine
   )

  parmFile = os.getcwd()+("/%sstack_"%prefixRun)+sample.fullName+".par"
  fp = file(parmFile, mode="w")
  fp.write(conf)

  # continue stacking if requested; pull files to local directory
  if os.access(treeFile, os.F_OK):
    os.system("mv %s %s" % (treeFile, "./"+tree.data))
    fp.write("loadSearchTree\n")
    fp.write("getTree\n")
    fp.write("doExtraction\n")
  else:
    fp.write("dumpSearchTree\n")
    fp.write("dumpTree\n")
    fp.write("doExtraction\n")
  fp.close()

  if not (continueRun and jobSuccessful(logFile, "Done!\n")):
    cmd = "%s --configFile=%s &> %s" % \
          (binPath, parmFile, logFile)
    os.system(cmd)

    if jobSuccessful(logFile, "Done!\n"):
      print "done"
    else:
      print "FAILED!"
      exit(-1)

  else:
    print "already done!"
    if os.access(parmFile, os.F_OK): os.unlink(parmFile)
    return

  minDist = None
  maxDist = None
  numVoids = None
  numParticles = None
  for line in open(logFile):
    if "Minimum void distance is" in line:
      line = line.split()
      minDist = float(line[4])/3000

    if "Maximum void distance is" in line:
      line = line.split()
      maxDist = float(line[4])/3000

    if "Selected" in line:
      line = line.split()
      numVoids = line[1]

    if "Number of particles =" in line:
      line = line.split()
      numParticles = line[4]

  open(voidDir+"/num_voids.txt", "w").write(numVoids)
  open(voidDir+"/num_particles.txt", "w").write(numParticles)

  if os.access(voidDir+"/NOVOID", os.F_OK):
    os.unlink(voidDir+"/NOVOID")

  if (numVoids == "0"):
    print "    No voids found; skipping!"
    fp = open(voidDir+"/NOVOID", "w")
    fp.write("no voids found\n")
    fp.close()
    return

  emptyStack = False
  for line in open(logFile):
    if "EMPTY STACK" in line:
      emptyStack = True
  if emptyStack:
    print "    Stack is empty; skipping!"
    fp = open(voidDir+"/NOVOID", "w")
    fp.write("empty stack\n")
    fp.close()
    return

  # figure out box volume and average density
  if sample.dataType == "observation":
    maskFile = sample.maskFile
    sulFunFile = sample.selFunFile

    if not os.access(sample.selFunFile, os.F_OK) and not sample.volumeLimited:
      print " Cannot find", selFunFile, "!"
      exit(-1)

    sys.stdout = open(logFile, 'a')
    sys.stderr = open(logFile, 'a')
    zMin = sample.zRange[0]
    zMax = sample.zRange[1]
    if not sample.volumeLimited:
      props = vp.getSurveyProps(maskFile, stack.zMin,
                                stack.zMax, zMin, zMax, "all",
                                selectionFuncFile=sample.selFunFile)
    else:
      zMinForVol = sample.zBoundary[0]
      zMaxForVol = sample.zBoundary[1]
      props = vp.getSurveyProps(maskFile, zMinForVol,
                                zMaxForVol, zMin, zMax, "all")
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
  
    boxVol = props[0]
    nbar   = props[1]
    if sample.volumeLimited:
      nbar = 1.0
  else:
    nbar = 1.0
 
    iX = float(sample.mySubvolume[0])
    iY = float(sample.mySubvolume[1])
    xMin = iX/sample.numSubvolumes * sample.boxLen
    yMin = iY/sample.numSubvolumes * sample.boxLen
    xMax = (iX+1)/sample.numSubvolumes * sample.boxLen
    yMax = (iY+1)/sample.numSubvolumes * sample.boxLen
    zMin = sample.zBoundaryMpc[0]
    zMax = sample.zBoundaryMpc[1]
    
    boxVol = (xMax-xMin)*(yMax-yMin)*(zMax-zMin)

  summaryLine = runSuffix + " " + \
                thisDataPortion + " " + \
                str(stack.zMin) + " " + \
                str(stack.zMax) + \
                " " + str(numVoids) + " " + str(minDist) + \
                " " + str(maxDist)
  if summaryFile != None:
    open(summaryFile, "a").write(summaryLine+"\n")

  # move local outputs to storage directory
  if os.access("tree.data", os.F_OK):
    normalization = float(numParticles)/float(boxVol)/nbar
    dataTemp = file("centers.txt", "r").readlines()
    fp = file("normalizations.txt", "w")
    for iVoid in xrange(len(dataTemp)):
      fp.write(str(normalization)+"\n")
    fp.close()

    os.system("mv %s %s" % ("tree.data", treeFile))
    os.system("mv %s %s" % ("void_indexes.txt", voidDir+"/"))
    os.system("mv %s %s" % ("posx.nc", voidDir+"/"))
    os.system("mv %s %s" % ("posy.nc", voidDir+"/"))
    os.system("mv %s %s" % ("posz.nc", voidDir+"/"))
    os.system("mv %s %s" % ("z_void_indexes.txt", voidDir+"/"))
    os.system("mv %s %s" % ("z_posx.nc", voidDir+"/"))
    os.system("mv %s %s" % ("z_posy.nc", voidDir+"/"))
    os.system("mv %s %s" % ("z_posz.nc", voidDir+"/"))
    os.system("mv %s %s" % ("redshifts.nc", voidDir+"/"))
    os.system("mv %s %s" % ("indexes.nc", voidDir+"/"))
    os.system("mv %s %s" % ("kdtree_stackvoids.dat", voidDir+"/"))
    os.system("mv %s %s" % ("centers.txt", voidDir+"/"))
    os.system("mv %s %s" % ("z_centers.txt", voidDir+"/"))
    os.system("mv %s %s" % ("sky_positions.txt", voidDir+"/"))
    os.system("mv %s %s" % ("check.txt", voidDir+"/"))
    os.system("mv %s %s" % ("tracer.txt", voidDir+"/"))
    os.system("mv %s %s" % ("normalizations.txt", voidDir+"/"))
    os.system("mv %s %s" % ("boundaryDistances.txt", voidDir+"/"))

  if os.access(idListFile, os.F_OK): os.unlink(idListFile)

  if os.access(parmFile, os.F_OK): os.unlink(parmFile)

  return
 
# -----------------------------------------------------------------------------
def launchCombine(sample, stack, voidDir=None, logFile=None,
                  zobovDir=None, workDir=None, thisDataPortion=None):

    sampleName = sample.fullName

    runSuffix = getStackSuffix(stack.zMin, stack.zMax, stack.rMin,
                               stack.rMax, thisDataPortion)

    sys.stdout = open(logFile, 'w')
    sys.stderr = open(logFile, 'a')

    if os.access(voidDir+"/num_voids.txt", os.F_OK):
      os.unlink(voidDir+"/num_voids.txt")

    doneGalUpdate = os.access(zobovDir+"donegalupdate", os.F_OK)

    for comboName in sample.comboList:
      if not os.access(zobovDir+"/galaxies.txt", os.F_OK):
        shutil.copy(workDir+"/sample_"+comboName+"/galaxies.txt", zobovDir)
      elif not doneGalUpdate:
        dataTemp = file(workDir+"/sample_"+comboName+"/galaxies.txt", 
                        "r").read()
        file(zobovDir+"/galaxies.txt", "a").write(dataTemp)

      sourceStackDir = workDir+"/sample_"+comboName+"/stacks_"+\
                       runSuffix

      if os.access(sourceStackDir+"/NOVOID", os.F_OK):
        continue

      if not os.access(voidDir+"/num_voids.txt", os.F_OK):
        # just copy
        shutil.copy(sourceStackDir+"/num_voids.txt", voidDir)
        shutil.copy(sourceStackDir+"/num_particles.txt", voidDir)
        shutil.copy(sourceStackDir+"/posx.nc", voidDir)
        shutil.copy(sourceStackDir+"/posy.nc", voidDir)
        shutil.copy(sourceStackDir+"/posz.nc", voidDir)
        shutil.copy(sourceStackDir+"/z_posx.nc", voidDir)
        shutil.copy(sourceStackDir+"/z_posy.nc", voidDir)
        shutil.copy(sourceStackDir+"/z_posz.nc", voidDir)
        shutil.copy(sourceStackDir+"/indexes.nc", voidDir)
        shutil.copy(sourceStackDir+"/redshifts.nc", voidDir)
        shutil.copy(sourceStackDir+"/centers.txt", voidDir)
        shutil.copy(sourceStackDir+"/void_indexes.txt", voidDir)
        shutil.copy(sourceStackDir+"/z_centers.txt", voidDir)
        shutil.copy(sourceStackDir+"/z_void_indexes.txt", voidDir)
        shutil.copy(sourceStackDir+"/sky_positions.txt", voidDir)
        shutil.copy(sourceStackDir+"/normalizations.txt", voidDir)
        shutil.copy(sourceStackDir+"/boundaryDistances.txt", voidDir)

      else:
        # append data
        dataTemp = file(voidDir+"/num_voids.txt", "r").readlines()
        dataTemp2 = file(sourceStackDir+"/num_voids.txt", "r").\
                    readlines()
        dataTemp  = np.array(dataTemp, dtype='i')
        dataTemp2 = np.array(dataTemp2, dtype='i')
        dataTemp += dataTemp2
        dataTemp = str(dataTemp[0])
        file(voidDir+"/num_voids.txt", "w").write(str(dataTemp))

        dataTemp = file(voidDir+"/num_particles.txt", "r").readlines()
        dataTemp2 = file(sourceStackDir+"/num_particles.txt", "r").\
                    readlines()
        dataTemp  = np.array(dataTemp, dtype='i')
        dataTemp2 = np.array(dataTemp2, dtype='i')
        dataTemp += dataTemp2
        dataTemp = str(dataTemp[0])
        file(voidDir+"/num_particles.txt", "w").write(str(dataTemp))

        dataTemp = file(sourceStackDir+"/z_centers.txt", "r").read()
        file(voidDir+"/z_centers.txt", "a").write(dataTemp)

        dataTemp = file(sourceStackDir+"/centers.txt", "r").read()
        file(voidDir+"/centers.txt", "a").write(dataTemp)

        dataTemp = file(sourceStackDir+"/normalizations.txt", "r").read()
        file(voidDir+"/normalizations.txt", "a").write(dataTemp)

        dataTemp = file(sourceStackDir+"/boundaryDistances.txt","r").read()
        file(voidDir+"/boundaryDistances.txt", "a").write(dataTemp)

        dataTemp = file(sourceStackDir+"/sky_positions.txt", "r").\
                   read()
        file(voidDir+"/sky_positions.txt", "a").write(dataTemp)

        idxTemp = file(sourceStackDir+"/z_void_indexes.txt", "r").\
                  readlines()
        idxTemp = np.array(idxTemp, dtype='i')
        fp = NetCDFFile(voidDir+"/z_posx.nc")
        dataTemp = fp.variables['array'][0:]
        fp.close()
        idxTemp[:] += len(dataTemp)
        fp = open(voidDir+"/z_void_indexes.txt", "a")
        for idx in idxTemp:
          fp.write(str(idx)+"\n")
        fp.close()
        dataTemp = ()

        fp = NetCDFFile(voidDir+"/z_posx.nc")
        dataTemp = fp.variables['array'][0:]
        fp.close()
        fp = NetCDFFile(sourceStackDir+"/z_posx.nc")
        dataTemp2 = fp.variables['array'][0:]
        fp.close()
        dataTemp = np.append(dataTemp, dataTemp2)
        outFile = NetCDFFile(voidDir+"/z_posx.nc", mode='w')
        outFile.createDimension("dim", len(dataTemp))
        v = outFile.createVariable("array", ncFloat, ("dim",))
        v[:] = dataTemp
        outFile.close()

        fp = NetCDFFile(voidDir+"/z_posy.nc")
        dataTemp = fp.variables['array'][0:]
        fp.close()
        fp = NetCDFFile(sourceStackDir+"/z_posy.nc")
        dataTemp2 = fp.variables['array'][0:]
        fp.close()
        dataTemp = np.append(dataTemp, dataTemp2)
        outFile = NetCDFFile(voidDir+"/z_posy.nc", mode='w')
        outFile.createDimension("dim", len(dataTemp))
        v = outFile.createVariable("array", ncFloat, ("dim",))
        v[:] = dataTemp
        outFile.close()

        fp = NetCDFFile(voidDir+"/z_posz.nc")
        dataTemp = fp.variables['array'][0:]
        fp.close()
        fp = NetCDFFile(sourceStackDir+"/z_posz.nc")
        dataTemp2 = fp.variables['array'][0:]
        fp.close()
        dataTemp = np.append(dataTemp, dataTemp2)
        outFile = NetCDFFile(voidDir+"/z_posz.nc", mode='w')
        outFile.createDimension("dim", len(dataTemp))
        v = outFile.createVariable("array", ncFloat, ("dim",))
        v[:] = dataTemp
        outFile.close()

        idxTemp = file(sourceStackDir+"/void_indexes.txt", "r").\
                  readlines()
        idxTemp = np.array(idxTemp, dtype='i')
        fp = NetCDFFile(voidDir+"/posx.nc")
        dataTemp = fp.variables['array'][0:]
        fp.close()
        idxTemp[:] += len(dataTemp)
        fp = open(voidDir+"/void_indexes.txt", "a")
        for idx in idxTemp:
          fp.write(str(idx)+"\n")
        fp.close()
        dataTemp = ()

        fp = NetCDFFile(voidDir+"/posx.nc")
        dataTemp = fp.variables['array'][0:]
        fp.close()
        fp = NetCDFFile(sourceStackDir+"/posx.nc")
        dataTemp2 = fp.variables['array'][0:]
        fp.close()
        dataTemp = np.append(dataTemp, dataTemp2)
        outFile = NetCDFFile(voidDir+"/posx.nc", mode='w')
        outFile.createDimension("dim", len(dataTemp))
        v = outFile.createVariable("array", ncFloat, ("dim",))
        v[:] = dataTemp
        outFile.close()

        fp = NetCDFFile(voidDir+"/posy.nc")
        dataTemp = fp.variables['array'][0:]
        fp.close()
        fp = NetCDFFile(sourceStackDir+"/posy.nc")
        dataTemp2 = fp.variables['array'][0:]
        fp.close()
        dataTemp = np.append(dataTemp, dataTemp2)
        outFile = NetCDFFile(voidDir+"/posy.nc", mode='w')
        outFile.createDimension("dim", len(dataTemp))
        v = outFile.createVariable("array", ncFloat, ("dim",))
        v[:] = dataTemp
        outFile.close()

        fp = NetCDFFile(voidDir+"/posz.nc")
        dataTemp = fp.variables['array'][0:]
        fp.close()
        fp = NetCDFFile(sourceStackDir+"/posz.nc")
        dataTemp2 = fp.variables['array'][0:]
        fp.close()
        dataTemp = np.append(dataTemp, dataTemp2)
        outFile = NetCDFFile(voidDir+"/posz.nc", mode='w')
        outFile.createDimension("dim", len(dataTemp))
        v = outFile.createVariable("array", ncFloat, ("dim",))
        v[:] = dataTemp
        outFile.close()


        fp = NetCDFFile(voidDir+"/redshifts.nc")
        dataTemp = fp.variables['array'][0:]
        fp.close()
        fp = NetCDFFile(sourceStackDir+"/redshifts.nc")
        dataTemp2 = fp.variables['array'][0:]
        fp.close()
        dataTemp = np.append(dataTemp, dataTemp2)
        outFile = NetCDFFile(voidDir+"/redshifts.nc", mode='w')
        outFile.createDimension("dim", len(dataTemp))
        v = outFile.createVariable("array", ncFloat, ("dim",))
        v[:] = dataTemp
        outFile.close()

        fp = NetCDFFile(voidDir+"/indexes.nc")
        dataTemp = fp.variables['array'][0:]
        fp.close()
        fp = NetCDFFile(sourceStackDir+"/indexes.nc")
        dataTemp2 = fp.variables['array'][0:]
        fp.close()
        dataTemp = np.append(dataTemp, dataTemp2)
        outFile = NetCDFFile(voidDir+"/indexes.nc", mode='w')
        outFile.createDimension("dim", len(dataTemp))
        v = outFile.createVariable("array", ncFloat, ("dim",))
        v[:] = dataTemp
        outFile.close()

    # add NOVOID if necessary
    if not os.access(voidDir+"/centers.txt", os.F_OK):
      fp = open(voidDir+"/NOVOID", "w")
      fp.write("no voids in combo\n")
      fp.close()

    fp = open(zobovDir+"/donegalupdate", "w")
    fp.close()

    print "Done!"

    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__

    if jobSuccessful(logFile, "Done!\n"):
      print "done"
    else:
      print "FAILED!"
      exit(-1)

  #else:
  #  print "already done!"

# -----------------------------------------------------------------------------
def launchProfile(sample, stack, voidDir=None, logFile=None, continueRun=None):

  sampleName = sample.fullName

  if os.access(voidDir+"/NOVOID", os.F_OK):
    print "no stack here; skipping!"
    return

  numVoids = open(voidDir+"/num_voids.txt", "r").readline()
  numParticles = open(voidDir+"/num_particles.txt", "r").readline()

  Rextracut = stack.rMin*3 + 1
  Rcircular = stack.rMin*3 + 2

  if not (continueRun and jobSuccessful(logFile, "Done!\n")):
    profileParmFile = voidDir+"/params.txt"
    fp = open(profileParmFile, "w")
    fp.write(str(stack.rMin)+"\n")
    fp.write(str(stack.rMax)+"\n")
    fp.write(str(Rcircular)+"\n")
    fp.write(numVoids+"\n")
    fp.close()

    #if sample.zRange[0] > stack.zMax or sample.zRange[1] < stack.zMin:
    #  print "outside sample redshift range; skipping!"
    #  fp = open(voidDir+"/NOVOID", "w")
    #  fp.write("outside z range: profile\n")
    #  fp.close()
    #  return

    if sample.profileBinSize == "auto":
      density = 0.5 * 50 / Rcircular / 2
    else:
      density = sample.profileBinSize

    sys.stdout = open(logFile, 'w')
    sys.stderr = open(logFile, 'a')
    vp.build_1d_profile(base_dir=voidDir, density=density,
                        rescaleMode=stack.rescaleMode)
    vp.build_2d_profile(base_dir=voidDir, density=density)
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__

    if jobSuccessful(logFile, "Done!\n"):
      print "done", numVoids
    else:
      print "FAILED!"
      exit(-1)

  else:
    print "already done!"


# -----------------------------------------------------------------------------
def launchFit(sample, stack, logFile=None, voidDir=None, figDir=None,
              continueRun=None, thisDataPortion=None):

  sampleName = sample.fullName

  runSuffix = getStackSuffix(stack.zMin, stack.zMax, stack.rMin,
                             stack.rMax, thisDataPortion)

  if not (continueRun and jobSuccessful(logFile, "Done!\n")):
    if os.access(voidDir+"/NOVOID", os.F_OK):
      print "no voids here; skipping!"
      return

    numVoids = int(open(voidDir+"/num_voids.txt", "r").readline())
    if numVoids < 10:
      print "not enough voids to fit; skipping!"
      fp = open(voidDir+"/NOFIT", "w")
      fp.write("not enough voids: %d \n" % numVoids)
      fp.close()
      return

    #if stack.zMin < sample.zRange[0] or stack.zMax > sample.zRange[1]:
    #  print "outside sample redshift range; skipping!"
    #  return

    if sample.partOfCombo or not sample.includeInHubble:
      print "sample not needed for further analysis; skipping!"
      fp = open(voidDir+"/NOFIT", "w")
      fp.write("sample not needed for hubble\n")
      fp.close()
      return

    if not stack.includeInHubble:
      print "radius not needed for further analysis; skipping!"
      return

    if not stack.includeInHubble:
      print "redshift not needed for further analysis; skipping!"
      return

    if os.access(figDir+"/stackedVoid_"+sampleName+"_"+\
                 runSuffix+".eps", os.F_OK):
      os.unlink(figDir+"/stackedVoid_"+sampleName+"_"+\
                runSuffix+".eps")

    sys.stdout = open(logFile, 'w')
    sys.stderr = open(logFile, 'a')

    badChain = True
    ntries = 0
    maxtries = 5
    while badChain:
      Rexpect = (stack.rMin+stack.rMax)/2
      Rtruncate = stack.rMin*3. + 1 # TEST
      #if sample.dataType == "observation":
      #  ret,fits,args = vp.fit_ellipticity(voidDir,Rbase=Rexpect,
      #                                Niter=300000,
      #                                Nburn=100000,
      #                                Rextracut=Rtruncate)
      #else:
      #  ret,fits,args = vp.fit_ellipticity(voidDir,Rbase=Rexpect,
      #                                Niter=300000,
      #                                Nburn=100000,
      #                                Rextracut=Rtruncate)
      #badChain = (args[0][0] > 0.5 or args[0][1] > stack.rMax or \
      #            args[0][2] > stack.rMax) and \
      #           (ntries < maxtries)
      ret,fits,args = vp.compute_inertia(voidDir, stack.rMax, mode="symmetric", nBootstraps=500)
      badChain = False
      ntries += 1

    #np.save(voidDir+"/chain.npy", ret)
    np.savetxt(voidDir+"/fits.out", fits)

    plotTitle = "Sample: "+sample.nickName+\
                ", z = "+str(stack.zMin)+"-"+\
                str(stack.zMax)+", R = "+str(stack.rMin)+"-"+\
                str(stack.rMax)+r" $h^{-1}$ Mpc"

    showRatio    = (ntries <= maxtries)
    showContours = (ntries <= maxtries)

    if stack.rescaleMode == "rmax":
      rescaleFactor = stack.rMax
    else:
      rescaleFactor = 1.0

    vp.draw_fit(*args, delta_tick=0.2, vmax=1.0, delta_v=0.01,
                nocontour=True, model_max=1.0, delta_model=0.1,
                plotTitle=plotTitle, show_ratio=showRatio,
                showContours=showContours,
                rescaleFactor=rescaleFactor)
    figure(1).savefig(figDir+"/stackedVoid_"+sampleName+\
                      "_"+runSuffix+".eps")
    figure(1).savefig(figDir+"/stackedVoid_"+sampleName+\
                      "_"+runSuffix+".pdf")
    figure(1).savefig(figDir+"/stackedVoid_"+sampleName+\
                      "_"+runSuffix+".png")

    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    if jobSuccessful(logFile, "Done!\n"):
      print "done (", ntries, " tries)"
    else:
      print "FAILED!"
      exit(-1)

    # record the measured stretch  
    exp = args[1]
    expError = args[2]
    outline = str(exp) + " " + str(expError)
    open(voidDir+"/expansion.out", "w").write(outline)

    if os.access(voidDir+"/NOFIT", os.F_OK):
      os.unlink(voidDir+"/NOFIT")
    if ntries > maxtries:
      print "    No reliable fit found; skipping!"
      fp = open(voidDir+"/NOFIT", "w")
      fp.write("bad ellipticity fit\n")
      fp.close()
      return

  else:
    print "already done!"


# -----------------------------------------------------------------------------
def launchHubble(dataPortions=None, dataSampleList=None, logDir=None,
                 INCOHERENT=None, workDir=None, figDir=None, errorBars=None, 
                 ZOBOV_PATH=None, continueRun=None, voidDir=None, 
                 doPlot = True, setName=None):

  for thisDataPortion in dataPortions:
    print "    For data portion", thisDataPortion
    sys.stdout.flush()

    # merge the bins from all samples
    numSamplesForHubble = 0
    maxRBins = 0
    maxZBins = 0
    allRBins = []
    allZBins = []

    for sample in dataSampleList:
      if sample.includeInHubble:
        numSamplesForHubble += 1

        for stack in sample.getHubbleStacks():
          alreadyHere = False
          for stackCheck in allRBins:
            if stack.rMin == stackCheck.rMin and stack.rMax==stackCheck.rMax:
              alreadyHere = True
              break
          if not alreadyHere:
            allRBins.append(stack)

          alreadyHere = False
          for stackCheck in allZBins:
            if stack.zMin == stackCheck.zMin and stack.zMax==stackCheck.zMax:
              alreadyHere = True
              break
          if not alreadyHere:
            allZBins.append(stack)

    allExpList = np.zeros((numSamplesForHubble,len(allRBins),len(allZBins),3))
    allExpList[:] = -1
    aveDistList = np.zeros((len(allZBins),3))

    iSample = -1
    for sample in dataSampleList:
      if not sample.includeInHubble:
        continue

      iSample += 1

      zobovDir = sample.zobovDir
      sampleName = sample.fullName

      print "      For data set", sampleName, "...",
      sys.stdout.flush()

      logFile = logDir+"/hubble_"+sampleName+"_"+thisDataPortion+".out"

      #AVEDIST_PATH = ZOBOV_PATH+"/mytools/computeAverageDistortionPMS"

      # compute the average stretch in each redshift bin for cleaner-looking
      #   (not not fully accurate) plots
      for (iZBin, stack) in enumerate(sample.getUniqueZBins()):

        if sample.dataType == "observation":
          aveDist = vp.aveStretchCone(stack.zMin, stack.zMax, 
                                      skyFrac = sample.skyFraction)
        else:
          aveDist = vp.aveStretch(stack.zMin, stack.zMax)

        aveDistList[iZBin, 0] = stack.zMin
        aveDistList[iZBin, 1] = aveDist
        aveDistList[iZBin, 2] = 0.00125

      numFits = 0
      fp = file(zobovDir+'/calculatedExpansions_'+thisDataPortion+'.txt',
                mode="w")

      expList = np.zeros((1, len(sample.getUniqueRBins()),
                         len(sample.getUniqueZBins()), 3))

      for (iZBin, zBin) in enumerate(sample.getUniqueZBins()):
        for (iR, rBin) in enumerate(sample.getUniqueRBins()):

          runSuffix = getStackSuffix(zBin.zMin, zBin.zMax, rBin.rMin,
                                     rBin.rMax, thisDataPortion)

          expFile = zobovDir+"/stacks_"+runSuffix+"/expansion.out"
          if not (os.access(zobovDir+"/stacks_"+runSuffix+"/NOVOID", \
                  os.F_OK) or
                  os.access(zobovDir+"/stacks_"+runSuffix+"/NOFIT", \
                  os.F_OK) or
              not os.access(expFile, os.F_OK)):
            exp      = float(open(expFile, "r").readline().split()[0])
            expError = float(open(expFile, "r").readline().split()[1])
          else:
            exp = -1
            expError = 0

          expList[0, iR, iZBin, 0] = exp
          expList[0, iR, iZBin, 1] = expError

          # TODO
          # compute the per-stack stretch for liklihood calculation
          #runSuffix = getStackSuffix(stack.zMin, stack.zMax, stack.rMin,
          #                            stack.rMax, thisDataPortion)
          #voidDir = zobovDir+"/stacks_" + runSuffix
          #centersFile = voidDir+"/centers.txt"
          #voidRedshifts = np.loadtxt(centersFile)[:,5]
          #aveDist = vp.aveWeightedStretchCone(stack.zMin, stack.zMax, 
          #                            skyFrac = sample.skyFraction, 
          #                            voidRedshifts=voidRedshifts)

          if sample.dataType == "observation":
            aveDist = vp.aveStretchCone(zBin.zMin, zBin.zMax, 
                                        skyFrac = sample.skyFraction)
          else:
            aveDist = vp.aveStretch(zBin.zMin, zBin.zMax)

          expList[0, iR, iZBin, 2] = aveDist

          for (iZCheck,zBinCheck) in enumerate(allZBins):
            for (iRCheck,rBinCheck) in enumerate(allRBins):
              if zBin.zMin == zBinCheck.zMin and zBin.zMax == zBinCheck.zMax:
                if rBin.rMin == rBinCheck.rMin and rBin.rMax == rBinCheck.rMax:
                  if sample.dataType == "observation":
                    aveDist = vp.aveStretchCone(zBin.zMin, zBin.zMax, 
                                                skyFrac = sample.skyFraction)
                  else:
                    aveDist = vp.aveStretch(zBin.zMin, zBin.zMax)

                  allExpList[iSample, iRCheck, iZCheck, 0] = exp
                  allExpList[iSample, iRCheck, iZCheck, 1] = expError

          fp.write(str(exp) + " " + str(expError) + " "+ " "+ \
                   str(aveDist) + " -1 " + runSuffix+"\n")
          numFits += 1

      fp.close()

      rlist = np.zeros((len(sample.getUniqueRBins()),2))
      for (iR,rBin) in enumerate(sample.getUniqueRBins()):
        rlist[iR, 0] = rBin.rMin
        rlist[iR, 1] = rBin.rMax

      zbase = np.zeros((len(sample.getUniqueZBins()),2))
      for (iZ,zBin) in enumerate(sample.getUniqueZBins()):
        zbase[iZBin,0] = zBin.zMinPlot
        zbase[iZBin,1] = zBin.zMaxPlot
      #np.savetxt(zobovDir+'/zbase_'+thisDataPortion+'.txt', zbase)

      if not (continueRun and jobSuccessful(logFile, "Done!\n")):
        sys.stdout = open(logFile, 'w')
        sys.stderr = open(logFile, 'a')
        plotTitle = "Sample: "+sample.nickName+", "+thisDataPortion+" voids"
        if doPlot:
          #vp.do_all_obs(zbase, expList, workDir+"/avedistortion_",
          vp.do_all_obs(zbase, expList, aveDistList,
                        rlist, plotTitle=plotTitle, plotAve=True)
          figure(1).savefig(figDir+"/hubble_"+sampleName+"_"+thisDataPortion+\
                            ".eps",bbox_inches='tight')
          figure(1).savefig(figDir+"/hubble_"+sampleName+"_"+thisDataPortion+\
                            ".pdf",bbox_inches='tight')
          figure(1).savefig(figDir+"/hubble_"+sampleName+"_"+thisDataPortion+\
                            ".png",bbox_inches='tight')
        else:
          print "Skipping plot"
          print "Done!"
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

        if jobSuccessful(logFile, "Done!\n"):
          print "done"
        else:
          print "FAILED!"
          exit(-1)

      else:
        print "already done!"

   # now do all samples in one plot 
    print "      For data set combined...",
    sys.stdout.flush()

    logFile = logDir+"/hubble_combined_"+thisDataPortion+".out"

    if not (continueRun and jobSuccessful(logFile, "Done!\n")):

      if errorBars == "ESTIMATED":
        # replace calculated errors with estimated ones
        errorBarFile = workDir+"/calculatedErrorBars_"+thisDataPortion+".txt"
        if os.access(errorBarFile, os.F_OK):
          print "REPLACING ERROR BARS!",
          sys.stdout.flush()
          fp = file(errorBarFile, mode="r")
          ignoreLine = fp.readline()
          for iZ in xrange(len(allExpList[0,0,:,0])):
            ignoreLine = fp.readline()

          for iZ in xrange(len(allExpList[0,0,:,0])):
            for iSample in xrange(len(allExpList[:,0,0,0])):
              for iR in xrange(len(allExpList[0,:,0,0])):
                line = fp.readline().split()
                allExpList[iSample,iR,iZ,2] = allExpList[iSample,iR,iZ,1]
                allExpList[iSample,iR,iZ,1] = float(line[1])

          fp.close()

      rlist = np.zeros((len(allRBins),2))
      for (iR,rBin) in enumerate(allRBins):
        rlist[iR, 0] = rBin.rMin
        rlist[iR, 1] = rBin.rMax

      zbase = np.zeros((len(allZBins),2))
      plotZmax = 0.0
      plotZmin = 1.e99
      for (iZ,zBin) in enumerate(allZBins):
        zbase[iZ,0] = zBin.zMinPlot
        zbase[iZ,1] = zBin.zMaxPlot
        
        if zBin.zMaxPlot > plotZmax: plotZmax = zBin.zMaxPlot
        if zBin.zMinPlot < plotZmin: plotZmin = zBin.zMinPlot

        if sample.dataType == "observation":
          aveDist = vp.aveStretchCone(zBin.zMin, zBin.zMax, 
                                      skyFrac = sample.skyFraction)
        else:
          aveDist = vp.aveStretch(zBin.zMin, zBin.zMax)

        aveDistList[iZ, 0] = zBin.zMin
        aveDistList[iZ, 1] = aveDist
        aveDistList[iZ, 2] = 0.00125


      shortSampleNames = list()
      for sample in dataSampleList:
        if sample.includeInHubble:
          shortSampleNames.append(sample.nickName)
      sys.stdout = open(logFile, 'w')
      sys.stderr = open(logFile, 'a')
      if doPlot:
        print "DOING PLOT"
        if INCOHERENT:
          #plotTitle = "all samples, incoherent "+\
          #            thisDataPortion+" voids"
          plotTitle = ''
        else:
          #plotTitle = "all samples, "+thisDataPortion+" voids"
          plotTitle = setName
          #plotTitle = ''
        vp.do_all_obs(zbase, allExpList, aveDistList,
                      rlist, plotTitle=plotTitle, sampleNames=shortSampleNames,
                      plotAve=True, mulfac = 1.0, biasLine = 1.16, 
                      plotZmin=plotZmin, plotZmax=plotZmax+0.2)
        figure(1).savefig(figDir+"/hubble_combined_"+thisDataPortion+\
                          ".eps",bbox_inches='tight')
        figure(1).savefig(figDir+"/hubble_combined_"+thisDataPortion+\
                          ".pdf",bbox_inches='tight')
        figure(1).savefig(figDir+"/hubble_combined_"+thisDataPortion+\
                          ".png",bbox_inches='tight')

        if INCOHERENT:
          plotTitle = "all samples, incoherent "+\
                      thisDataPortion+" voids (systematics corrected)"
        else:
          #plotTitle = "all samples, "+thisDataPortion+\
          #            " voids (systematics corrected)"
          plotTitle = setName + "(sysematics corrected)"
        vp.do_all_obs(zbase, allExpList, aveDistList,
                      rlist, plotTitle=plotTitle, sampleNames=shortSampleNames,
                      plotAve=True, mulfac = 1.16, 
                      plotZmin=plotZmin, plotZmax=plotZmax+0.2)
        figure(1).savefig(figDir+"/hubble_combined_"+thisDataPortion+\
                          "_debiased.eps",bbox_inches='tight')
        figure(1).savefig(figDir+"/hubble_combined_"+thisDataPortion+\
                          "_debiased.pdf",bbox_inches='tight')
        figure(1).savefig(figDir+"/hubble_combined_"+thisDataPortion+\
                          "_debiased.png",bbox_inches='tight')
      else:
        print "Skipping plot"
        print "Done!"
      sys.stdout = sys.__stdout__
      sys.stderr = sys.__stderr__

      # save all expansion data to a single file 
      fp = file(workDir+'/calculatedExpansions_'+thisDataPortion+'.txt',
                mode="w")
      fp.write(str(numSamplesForHubble) + " " +
               str(len(allRBins)) + " " + str(len(allZBins)) + "\n")
      for zBin in allZBins:
        fp.write(str(zBin.zMin) + " " + str(zBin.zMax) + "\n")

      for iZ in xrange(len(allExpList[0,0,:,0])):
        for iSample in xrange(len(allExpList[:,0,0,0])):
          for iR in xrange(len(allExpList[0,:,0,0])):
            fp.write(str(allExpList[iSample,iR,iZ,0]) + " " +\
                     str(allExpList[iSample,iR,iZ,1]) + " " +\
                     str(allExpList[iSample,iR,iZ,2]) + "\n")
      fp.close()

      # save all void distribution data to a single file 
      fp = file(workDir+'/voidDistributions_'+thisDataPortion+'.txt',
                mode="w")
      fp.write(str(numSamplesForHubble) + " " +
               str(len(allRBins)) + " " + str(len(allZBins)) + "\n")
      for zBin in allZBins:
        fp.write(str(zBin.zMin) + " " + str(zBin.zMax) + "\n")

      for zBin in allZBins:
        for sample in dataSampleList:
          if not sample.includeInHubble: continue
          for rBin in allRBins:
            runSuffix = getStackSuffix(zBin.zMin, zBin.zMax, rBin.rMin,
                                      rBin.rMax, thisDataPortion)
            voidDir = sample.zobovDir+"/stacks_" + runSuffix
            centersFile = voidDir+"/centers.txt"
            if os.access(centersFile, os.F_OK):
              voidRedshifts = np.loadtxt(centersFile)
              if voidRedshifts.ndim > 1:
                voidRedshifts = voidRedshifts[:,5]
              else:
                voidRedshifts = voidRedshifts[5]
              #fp.write(str(len(voidRedshifts))+" ")
              np.savetxt(fp, voidRedshifts[None])
            else:
              fp.write("-1\n")
      fp.close()


      if jobSuccessful(logFile, "Done!\n"):
        print "done"
      else:
        print "FAILED!"
        exit(-1)

    else:
      print "already done!"

# -----------------------------------------------------------------------------
def launchLikelihood(dataPortions=None, logDir=None, workDir=None, 
                      continueRun=False):

  for thisDataPortion in dataPortions:
    print "    For data portion", thisDataPortion, "...",
    sys.stdout.flush()

    logFile = logDir+"/likelihood_"+thisDataPortion+".out"
    
    if not (continueRun and jobSuccessful(logFile, "Done!\n")):

      sys.stdout = open(logFile, 'w')
      sys.stderr = open(logFile, 'a')

      vp.build1dLikelihood(workDir+"/calculatedExpansions_"+\
                         thisDataPortion+".txt",
                    workDir+"/voidDistributions_" + \
                      thisDataPortion+".txt",
                    nbins = 20,
                    OmStart = 0.0,
                    OmEnd   = 1.0,
                    biasStart = 1.0,
                    biasEnd   = 1.32,
                    #biasStart = 1.15,
                    #biasEnd   = 1.17,
                    outputBase = workDir+"/1dlikelihoods_"+thisDataPortion+"_",
                    useBinAve = False)

      sys.stdout = sys.__stdout__
      sys.stderr = sys.__stderr__

      if jobSuccessful(logFile, "Done!\n"):
        print "done"
      else:
        print "FAILED!"
        exit(-1)

    else:
      print "already done!"


# -----------------------------------------------------------------------------
def launchEstimateErrorBars(workDir=None, nullDir=None, numNulls=None,
                            dataPortion=None, nullDirList=None):

  IVALUE = 0
  IERROR = 1
  
  # open first null file to get sizes
  if nullDirList == None:
    nullFile = nullDir + "/" + str(1) + "/calculatedExpansions_" + \
               str(dataPortion) + ".txt"
  else:
    nullFile = nullDirList[0] + "/calculatedExpansions_" + str(dataPortion) + \
               ".txt"
    numNulls = len(nullDirList)

  fp = open(nullFile, "r")
  binLine = fp.readline()
  numSamples = int(binLine.split(" ")[0])
  numRBins = int(binLine.split(" ")[1])
  numZBins = int(binLine.split(" ")[2])
  numBins = numRBins * numSamples * numZBins
  nRTotal = numRBins * numSamples
  fp.close()
  
  # load null data
  nullData = np.zeros([numNulls, numBins, 2])
  
  for iNull in xrange(numNulls):
    if nullDirList == None:
      nullFile = nullDir + "/" + str(iNull) + "/calculatedExpansions_" + \
                 str(dataPortion) + ".txt"
    else:
      nullFile = nullDirList[iNull] + "/calculatedExpansions_" + \
                 str(dataPortion) + ".txt"

    fp = open(nullFile, "r")
  
    binLine = fp.readline()
    numSamples = int(binLine.split(" ")[0])
    numRBins = int(binLine.split(" ")[1])
    numZBins = int(binLine.split(" ")[2])
    numBins = numRBins * numSamples * numZBins
  
    zList = np.zeros( [numZBins,2] )
    for iZ in xrange(int(numZBins)):
      line = fp.readline().split()
      zList[iZ,0] = float(line[0])
      zList[iZ,1] = float(line[1])
  
    for iZ in xrange(numZBins):
      iRTotal = 0
      for iSample in xrange(numSamples):
        for iR in xrange(numRBins):
          line = fp.readline().split()
          nullData[iNull,iRTotal*numZBins+iZ,IVALUE] = float(line[0])
          nullData[iNull,iRTotal*numZBins+iZ,IERROR] = float(line[1])
          iRTotal += 1
    fp.close()
  
  # use 68% limits to get error bars
  errorBars = np.zeros([numBins, 2])
  for iBin in xrange(numBins):
    binData = nullData[:,iBin,IVALUE]
    validData = binData[:] > -1
    binData = binData[validData]
    binData = np.sort(binData)
  
    if np.size(binData) > 0:
      errorBars[iBin, 0] = binData[int(0.16*np.size(binData))]
      errorBars[iBin, 1] = binData[int(0.84*np.size(binData))]
 
      mean = (errorBars[iBin, 1] + errorBars[iBin, 0])/2
      sig =  (errorBars[iBin, 1] - errorBars[iBin, 0])/2
      errorBars[iBin, 0] = mean
      errorBars[iBin, 1] = sig
    else:
      errorBars[iBin, :] = -1
  
  # writes out error bars
  fp = file(workDir+"/calculatedErrorBars_" + \
            str(dataPortion) + ".txt", mode="w")
  fp.write(str(numSamples) + " " + str(numRBins) + " " + str(numZBins) + "\n")
  for iZ in xrange(numZBins):
    fp.write(str(zList[iZ,0]) + " " + str(zList[iZ,1]) + "\n")
  
  for iZ in xrange(numZBins):
    iRTotal = 0
    for iSample in xrange(numSamples):
      for iR in xrange(numRBins):
        fp.write(str(errorBars[iRTotal*numZBins+iZ,0]) + " " + \
                 str(errorBars[iRTotal*numZBins+iZ,1]) + "\n")
        iRTotal += 1
  fp.close()
  
  
                      

# -----------------------------------------------------------------------------
def launchVelocityStack(sample, stack, binPath, 
                         velField_file,
                         thisDataPortion=None, logDir=None,
                         voidDir=None, runSuffix=None,
                         zobovDir=None,
                         summaryFile=None, 
                         continueRun=None, dataType=None, prefixRun=""):

  sampleName = sample.fullName

  runSuffix = getStackSuffix(stack.zMin, stack.zMax, stack.rMin,
                             stack.rMax, thisDataPortion)

  logFile = logDir+("/%svelocity_stack_"%prefixRun)+sampleName+"_"+runSuffix+".out"

  voidCenters=voidDir+"/centers.txt"
#  Rmax = 

  centralRadius = stack.rMin * 0.25
  Rextracut = stack.rMin*3 + 1
  Rcircular = stack.rMin*3 + 2

  parameters="--velocityField=%s --voidCenters=%s --Rmax=%e --L0=%e --numBins=%d" % (velField_file, voidCenters, Rmax, Boxsize, numBins)

  if not (continueRun and jobSuccessful(logFile, "Done!\n")):
    cmd = "%s %s &> %s" % (binPath,parameters,logFile)
    os.system(cmd)
    if jobSuccessful(logFile, "Done!\n"):
      print "done"
    else:
      print "FAILED!"
      exit(-1)
      
  else:
    print "already done!"
