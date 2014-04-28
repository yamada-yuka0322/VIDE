#+
#   VIDE -- Void IDentification and Examination -- ./crossCompare/plotting/plotDenMaps.py
#   Copyright (C) 2010-2014 Guilhem Lavaux
#   Copyright (C) 2011-2014 P. M. Sutter
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
#!/usr/bin/env python
#+
#   VIDE -- Void IDentification and Examination -- ./pipeline/apAnalysis.py
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

# takes nVoids evenly distributed, plots a slice of the local density and
#  overlays the voids

import imp
import pickle
import os
import numpy as np
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from void_python_tools.backend import *
import void_python_tools.xcor as xcor
from netCDF4 import Dataset
#import pylab as plt

NetCDFFile = Dataset
ncFloat = 'f8'

matplotlib.rcParams.update({'font.size': 16})
# ------------------------------------------------------------------------------

mergerNameBase = "voidOverlap"

parser = argparse.ArgumentParser(description='Analyze.')
parser.add_argument('--parm', dest='parm', default='datasetsToAnalyze.py', help='path to parameter file')
parser.add_argument('--show', dest='showPlot', action='store_const',
                   const=True, default=False,
                   help='display the plot (default: just write eps)')
args = parser.parse_args()

nVoids = 10

# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# plot a slice of the density around the void in baseIDList, 
#   with any voids in the slice shown and any voids in baseIDList flagged
def plotVoidAndDen(idList, voidList, partData, boxLen, figDir,
                   sliceCenter=None, sliceWidth=200,
                   baseIDList=None, baseRadius=0, nickName=None,
                   baseNickName=None,
                   baseIndex=0,
                   periodic=None, showPlot=False, plotName=None):

  if len(voidList) <= 1: return
  plt.clf()

  #sliceWidth = 220
  sliceWidth = max(220, sliceWidth)

  # make an appropriate box
  xwidth = sliceWidth
  ywidth = sliceWidth
  zwidth = sliceWidth/4.
  #zwidth = max(sliceWidth/4., 50)

  # get mean density
  part = 1.*partData
  totalNumPart = len(part)
  totalVol = (part[:,0].max() - part[:,0].min()) * \
             (part[:,1].max() - part[:,1].min()) * \
             (part[:,2].max() - part[:,2].min())

  meanDen = totalNumPart/totalVol

  # single out the matched void
  keepVoid = []
  if len(np.atleast_1d(idList)) > 0:
    keepVoid = voidList[voidList[:,7] == idList]
  if len(np.shape(keepVoid)) > 1: keepVoid = keepVoid[0,:]
  filter = voidList[:,7] != idList
  voidList = voidList[filter,:]

  # convert everything to relative coordinates
  part[:,0] -= sliceCenter[0]
  part[:,1] -= sliceCenter[1]
  part[:,2] -= sliceCenter[2]

  shiftUs = np.abs(part[:,0]) > boxLen[0]/2.
  if ("x" in periodicLine): part[shiftUs,0] -= \
                            np.copysign(boxLen[0],part[shiftUs,0])
  shiftUs = np.abs(part[:,1]) > boxLen[1]/2.
  if ("y" in periodicLine): part[shiftUs,1] -= \
                            np.copysign(boxLen[1],part[shiftUs,1])
  shiftUs = np.abs(part[:,2]) > boxLen[2]/2.
  if ("z" in periodicLine): part[shiftUs,2] -= \
                            np.copysign(boxLen[2],part[shiftUs,2])
  
  voidList = np.atleast_2d(voidList)
  np.atleast_2d(voidList)[:,0] -= sliceCenter[0]
  np.atleast_2d(voidList)[:,1] -= sliceCenter[1]
  np.atleast_2d(voidList)[:,2] -= sliceCenter[2]

  shiftUs = np.abs(voidList[:,0]) > boxLen[0]/2.
  if ("x" in periodicLine): 
    voidList[shiftUs,0] -= \
                            np.copysign(boxLen[0],voidList[shiftUs,0])
  shiftUs = np.abs(voidList[:,1]) > boxLen[1]/2.
  if ("y" in periodicLine): voidList[shiftUs,1] -= \
                            np.copysign(boxLen[1],voidList[shiftUs,1])
  shiftUs = np.abs(voidList[:,2]) > boxLen[2]/2.
  if ("z" in periodicLine): voidList[shiftUs,2] -= \
                            np.copysign(boxLen[2],voidList[shiftUs,2])

  if len(np.atleast_1d(keepVoid)) >= 1:
    keepVoid[0] -= sliceCenter[0]
    keepVoid[1] -= sliceCenter[1]
    keepVoid[2] -= sliceCenter[2]
    shiftUs = np.abs(keepVoid[0]) > boxLen[0]/2.
    if ("x" in periodicLine) and shiftUs: keepVoid[0] -= \
                                        np.copysign(boxLen[0],keepVoid[0])
    shiftUs = np.abs(keepVoid[1]) > boxLen[1]/2.
    if ("y" in periodicLine) and shiftUs: keepVoid[1] -= \
                                        np.copysign(boxLen[1],keepVoid[1])
    shiftUs = np.abs(keepVoid[2]) > boxLen[2]/2.
    if ("z" in periodicLine) and shiftUs: keepVoid[2] -= \
                                        np.copysign(boxLen[2],keepVoid[2])

  xmin = -xwidth/2.
  xmax = xwidth/2.
  ymin = -ywidth/2.
  ymax = ywidth/2.
  zmin = -zwidth/2.
  zmax = zwidth/2.

  # pull out voids that were potential matches
  filter = np.sqrt(voidList[:,0]**2 + voidList[:,1]**2 + voidList[:,2]**2) <=\
           baseRadius*1.5
  potentialMatches = voidList[filter] 

  # get centers and radii of any other voids in slice
  zminVoid = -zwidth/16.
  zmaxVoid = zwidth/16.

  filter = (voidList[:,0] > xmin) & (voidList[:,0] < xmax) & \
           (voidList[:,1] > ymin) & (voidList[:,1] < ymax) & \
           (voidList[:,2] > zminVoid) & (voidList[:,2] < zmaxVoid) 
  voidList = voidList[filter,:]

  # slice particles
  filter = (part[:,0] > xmin) & (part[:,0] < xmax) & \
           (part[:,1] > ymin) & (part[:,1] < ymax) & \
           (part[:,2] > zmin) & (part[:,2] < zmax) 
  part = part[filter]

  # plot density
  extent = [xmin, xmax, ymin, ymax]
  hist, xedges, yedges = np.histogram2d(part[:,0], part[:,1], normed=False,
                                        bins=64)
 
  #hist /= meanDen 
  hist = np.log10(hist+1)
  plt.imshow(hist, 
             aspect='equal',
             extent=extent,
             interpolation='gaussian',
             cmap='YlGnBu_r')
  #plt.colorbar()

  # overlay voids as circles
  fig = plt.gcf()
  ax = fig.add_subplot(1,1,1)

  # the original void
  circle = plt.Circle((0,0), baseRadius,
                      edgecolor='orange', facecolor=None, fill=False,
                      linewidth=5)
  fig.gca().add_artist(circle)

  # our matched void
  if len(np.atleast_1d(keepVoid)) > 0:
    if idList == baseIDList:
      edgecolor = 'orange'
    else:
      edgecolor = 'red'

    circle = plt.Circle((keepVoid[0], keepVoid[1]), keepVoid[4], 
                         edgecolor=edgecolor, facecolor=None, fill=False,
                         linewidth=5)
    fig.gca().add_artist(circle)
 
  # other voids in the slice 
  for void in voidList:
    if np.any(void[7] == idList):
      continue
    else:
      color = 'white'
    circle = plt.Circle((void[0], void[1]), void[4], 
                         edgecolor=color, facecolor=None, fill=False,
                         linewidth=5)
    fig.gca().add_artist(circle)

  # potential match voids
  #for void in potentialMatches:
  #  if np.any(void[7] == idList):
  #    continue
  #  else:
  #    color = 'green'
  #  circle = plt.Circle((void[0], void[1]), void[4], 
  #                       edgecolor=color, facecolor=None, fill=False,
  #                       linewidth=5)
  #  fig.gca().add_artist(circle)

  baseNickName = baseNickName[:-10].lstrip()
  nickName = nickName[:-10].lstrip()
  if idList == baseIDList:
    title = "%d $h^{-1}$Mpc" % int(baseRadius)
    title += " (" + baseNickName + ")"
  else:
    title = r"$\rightarrow$ " 
    if len(np.atleast_1d(keepVoid)) > 0:
      title += "%d $h^{-1}$Mpc" % int(keepVoid[4]) + " (" + nickName + ")"
    else:
      title += "No match"
      
  #title += "(" + str(int(baseIDList)) + ")"
  
  plt.title(title, fontsize=20)
  #plt.xlabel("x [$h^{-1}$Mpc]", fontsize=14)
  #plt.ylabel("y [$h^{-1}$Mpc]", fontsize=14)

  plotName += "_" + str(int(baseIndex))
  #plotName += "_" + str(int(baseRadius))
  plt.savefig(figDir+"/fig_"+plotName+".pdf", bbox_inches="tight")
  plt.savefig(figDir+"/fig_"+plotName+".eps", bbox_inches="tight")
  plt.savefig(figDir+"/fig_"+plotName+".png", bbox_inches="tight")

  if showPlot: os.system("display %s" % figDir+"/fig_"+plotName+".png")

  return

# -----------------------------------------------------------------------------

filename = args.parm
print " Loading parameters from", filename
if not os.access(filename, os.F_OK):
  print "  Cannot find parameter file %s!" % filename
  exit(-1)
parms = imp.load_source("name", filename)
globals().update(vars(parms))

if not os.access(outputDir, os.F_OK):
  os.makedirs(outputDir)

if not os.access(logDir, os.F_OK):
  os.makedirs(logDir)

if not os.access(figDir, os.F_OK):
  os.makedirs(figDir)

mergerFileBase = outputDir + "/" + mergerNameBase

# get list of base voids
with open(workDir+baseSampleDir+"/sample_info.dat", 'rb') as input:
  baseSample = pickle.load(input)
  baseSampleName = baseSample.fullName
  baseVoidList = np.loadtxt(workDir+baseSampleDir+"/centers_central_"+\
                    baseSampleName+".out")

  # sort by size
  radii = baseVoidList[:,4]
  indices = np.argsort(radii)[::-1]
  baseVoidList = baseVoidList[indices,:]
  setName = baseSampleDir.split('/')[0]

# pick our void sample
bigVoidList = baseVoidList[0:10,:]
stride = len(baseVoidList)/nVoids
baseVoidList = baseVoidList[::stride]
baseVoidList = np.vstack((bigVoidList,baseVoidList))

sampleDirList.insert(0,baseSampleDir)

for (iSample, sampleDir) in enumerate(sampleDirList):
  if compareSampleTag in sampleDir: continue
  with open(workDir+sampleDir+"/sample_info.dat", 'rb') as input:
    sample = pickle.load(input)
  print "  Working with", sample.fullName, "..."
  sys.stdout.flush()
  sampleName = sample.fullName

  print "    Loading particle data..."
  sys.stdout.flush()

  
  infoFile = workDir+"/"+sampleDir+"/zobov_slice_"+sample.fullName+".par"
  File = NetCDFFile(infoFile, 'r')
  ranges = np.zeros((3,2))
  ranges[0][0] = getattr(File, 'range_x_min')
  ranges[0][1] = getattr(File, 'range_x_max')
  ranges[1][0] = getattr(File, 'range_y_min')
  ranges[1][1] = getattr(File, 'range_y_max')
  ranges[2][0] = getattr(File, 'range_z_min')
  ranges[2][1] = getattr(File, 'range_z_max')
  File.close()
  mul = np.zeros((3))
  mul[:] = ranges[:,1] - ranges[:,0]
  boxLen = mul

  partFile = workDir+"/"+sampleDir+"/zobov_slice_"+sample.fullName
  #partFile = catalogDir+"/"+sample.dataFile
  iLine = 0
  partData = []
  part = np.zeros((3))
  File = file(partFile)
  chk = np.fromfile(File, dtype=np.int32,count=1)
  Np = np.fromfile(File, dtype=np.int32,count=1)
  chk = np.fromfile(File, dtype=np.int32,count=1)

  chk = np.fromfile(File, dtype=np.int32,count=1)
  x = np.fromfile(File, dtype=np.float32,count=Np)
  x *= mul[0]
  x += ranges[0][0]
  chk = np.fromfile(File, dtype=np.int32,count=1)

  chk = np.fromfile(File, dtype=np.int32,count=1)
  y = np.fromfile(File, dtype=np.float32,count=Np)
  y *= mul[1] 
  y += ranges[1][0]
  chk = np.fromfile(File, dtype=np.int32,count=1)

  chk = np.fromfile(File, dtype=np.int32,count=1)
  z = np.fromfile(File, dtype=np.float32,count=Np)
  z *= mul[2] 
  z += ranges[2][0]
  chk = np.fromfile(File, dtype=np.int32,count=1)
  File.close()
  partData = np.column_stack((x,y,z))#.transpose()

  for (iBaseVoid,baseVoid) in enumerate(baseVoidList):
    print "   Void:", int(baseVoid[7]), "(", int(baseVoid[4]), ")"
    baseIDList = baseVoid[7]
    sliceCenter = baseVoid[0:3]
    sliceWidth  = baseVoid[4]*4
 
    # get matched void
    idList = []
    if sample.fullName == baseSample.fullName:
      idList = baseIDList
    else:
      matchFile=mergerFileBase+"_"+baseSampleName+"_"+sampleName+"_summary.out"
      if os.access(matchFile, os.F_OK):
        matchList = np.loadtxt(matchFile)
        for i,testID in enumerate(matchList[:,0]):
          if testID == baseIDList:
            if (matchList[i,8] > 0): idList.append(matchList[i,8])
        idList = np.array(idList)
        idList = idList.astype(int)

    voidList = np.loadtxt(workDir+sampleDir+"/trimmed_nodencut_centers_central_"+\
                          sampleName+".out")

    periodicLine = getPeriodic(sample)
    plotVoidAndDen(idList, voidList, partData, boxLen, figDir, 
                   sliceCenter=sliceCenter, sliceWidth=sliceWidth, 
                   baseIDList=baseIDList, baseRadius=baseVoid[4],
                   baseIndex=iBaseVoid,
                   nickName=sample.nickName, periodic=periodicLine,
                   baseNickName=baseSample.nickName,
                   showPlot=args.showPlot, 
                   plotName="denmap_"+setName+"_"+baseSampleName+"_"+sampleName) 
print " Done!"

