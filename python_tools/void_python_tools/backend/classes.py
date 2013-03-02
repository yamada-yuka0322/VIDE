#+
#   VIDE -- Void IDEntification pipeline -- ./python_tools/void_python_tools/backend/classes.py
#   Copyright (C) 2010-2013 Guilhem Lavaux
#   Copyright (C) 2011-2013 Paul M. Sutter
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
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
# classes and routines used to support scripts

import os

LIGHT_SPEED = 299792.458

class Stack:
  zMin = 0.0
  zMax = 0.1
  rMin = 5
  rMax = 15
  zMinPlot = 0.0
  zMaxPlot = 0.1
  includeInHubble = True
  partOfCombo = False
  needProfile = True
  rescaleMode = "rmax" # options: "rmax" to scale to largest void in stack
                       #          "rv"   normalize each void

  def __init__(self, zMin, zMax, rMin, rMax, includeInHubble, partOfCombo, 
               zMinPlot=None, needProfile=True, rescaleMode="rmax"):
    self.zMin = zMin
    self.zMax = zMax
    self.rMin = rMin
    self.rMax = rMax
    self.zMaxPlot = zMax
    self.includeInHubble = includeInHubble
    self.partOfCombo = partOfCombo
    self.needProfile = needProfile
    self.rescaleMode = rescaleMode
   
    if zMinPlot == None:
      self.zMinPlot = self.zMin
    else:
      self.zMinPlot = zMinPlot

class Sample:
  dataType = "observation"
  dataFormat = "sdss"
  dataFile = "lss.dr72dim.dat"
  dataUNit = 1
  fullName = "lss.dr72dim.dat"
  nickName = "dim"
  zobovDir = ""
  maskFile = "rast_window_512.fits"
  selFunFile = "czselfunc.all.dr72dim.dat"
  skyFraction = 0.19
  zBoundary = (0.0, 0.1)
  zBoundaryMpc = (0., 300)
  zRange    = (0.0, 0.1)
  omegaM    = 0.27
  minVoidRadius = 5
  fakeDensity = 0.01
  profileBinSize = 2 # Mpc
  volumeLimited = True
  includeInHubble = True
  partOfCombo  = False
  isCombo = False
  comboList = []

  # applies to simulations only
  boxLen = 1024 # Mpc/h
  usePecVel = False
  subsample = 1.0
  useLightCone = True
  numSubvolumes = 1
  mySubvolume = 1

  stacks = []

  def __init__(self, dataFile="", fullName="", dataUnit=1,
               nickName="", maskFile="", selFunFile="",
               zBoundary=(), zRange=(), zBoundaryMpc=(),
               minVoidRadius=0, fakeDensity=0.01, volumeLimited=True,
               includeInHubble=True, partOfCombo=False, isCombo=False, 
               comboList=(), profileBinSize=2.0, skyFraction=0.19,
               boxLen=1024, usePecVel=False, omegaM=0.27, 
               numSubvolumes=1, mySubvolume=1, dataFormat="sdss",
               dataType="observation",
               subsample=1.0, useLightCone=True):
    self.dataFile = dataFile
    self.fullName = fullName
    self.nickName = nickName
    self.maskFile = maskFile
    self.selFunFile = selFunFile
    self.zBoundary = zBoundary 
    self.zBoundaryMpc = zBoundaryMpc
    self.zRange = zRange
    self.minVoidRadius = minVoidRadius
    self.fakeDensity = fakeDensity
    self.volumeLimited = volumeLimited
    self.includeInHubble = includeInHubble 
    self.partOfCombo = partOfCombo 
    self.isCombo = isCombo
    self.comboList = comboList
    self.zobovDir = None
    self.profileBinSize = profileBinSize
    self.skyFraction = skyFraction
    self.dataType = dataType
    self.boxLen = boxLen
    self.usePecVel = usePecVel
    self.omegaM = omegaM
    self.numSubvolumes = numSubvolumes
    self.mySubvolume = mySubvolume
    self.dataFormat = dataFormat
    self.subsample = subsample
    self.useLightCone = useLightCone
    self.dataUnit = dataUnit

    self.stacks = []

  def addStack(self, zMin, zMax, rMin, rMax, 
                includeInHubble, partOfCombo,zMinPlot=None, 
                needProfile=True, rescaleMode="rmax"):
    self.stacks.append(Stack(zMin, zMax, rMin, rMax,
                            includeInHubble, partOfCombo, 
                            zMinPlot=zMinPlot, needProfile=needProfile,
                            rescaleMode=rescaleMode))

  def getHubbleStacks(self):
    stacksForHubble = []
    for stack in self.stacks:
      if stack.includeInHubble:
        stacksForHubble.append(stack)
    return stacksForHubble

  def getUniqueRBins(self):
    uniqueRStacks = []
    for stack in self.stacks:
      if not stack.includeInHubble:
        continue
      alreadyHere = False
      for stackCheck in uniqueRStacks:
        if stack.rMin == stackCheck.rMin and stack.rMax == stackCheck.rMax:
          alreadyHere = True
          break
      if not alreadyHere:
        uniqueRStacks.append(stack)
    return uniqueRStacks

  def getUniqueZBins(self):
    uniqueZStacks = []
    for stack in self.stacks:
      if not stack.includeInHubble:
        continue
      alreadyHere = False
      for stackCheck in uniqueZStacks:
        if stack.zMin == stackCheck.zMin and stack.zMax == stackCheck.zMax:
          alreadyHere = True
          break
      if not alreadyHere:
        uniqueZStacks.append(stack)
    return uniqueZStacks

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
def jobSuccessful(logFile, doneString):
  jobDone = False
  checkLine = ""
  if os.access(logFile, os.F_OK):
    #filelines = file(logFile, "r").readlines()
    #if len(filelines) >= 1:
    #  checkLine = filelines[-1]
    for line in open(logFile, 'r'):
      if doneString in line: jobDone = True
  return jobDone

def getStackSuffix(zMin, zMax, rMin, rMax, dataPortion):
  return "z"+str(zMin)+"-"+str(zMax)+"_"+str(rMin)+"-"+str(rMax)+\
         "Mpc"+"_"+dataPortion

def my_import(name):
    mod = __import__(name)
    components = name.split('.')
    for comp in components[1:]:
        mod = getattr(mod, comp)
    return mod
