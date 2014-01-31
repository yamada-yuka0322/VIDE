#+
#   VIDE -- Void IDentification and Examination -- ./python_tools/void_python_tools/plotting/__init__.py
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

# Various utilities for loading and modifying particle datasets

import numpy as np
from netCDF4 import Dataset
import sys
from void_python_tools.backend import *
import void_python_tools.apTools as vp

NetCDFFile = Dataset
ncFloat = 'f8'

# -----------------------------------------------------------------------------
def loadPart(workDir, sampleDir, sample):
    #print "    Loading particle data..."
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
    isObservation = getattr(File, 'is_observation')
    maskIndex = getattr(File, 'mask_index')
    File.close()
    mul = np.zeros((3))
    mul[:] = ranges[:,1] - ranges[:,0]

    partFile = workDir+"/"+sampleDir+"/zobov_slice_"+sample.fullName
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
    if isObservation != 1:
      x += ranges[0][0]
    chk = np.fromfile(File, dtype=np.int32,count=1)

    chk = np.fromfile(File, dtype=np.int32,count=1)
    y = np.fromfile(File, dtype=np.float32,count=Np)
    y *= mul[1] 
    if isObservation != 1:
      y += ranges[1][0]
    chk = np.fromfile(File, dtype=np.int32,count=1)

    chk = np.fromfile(File, dtype=np.int32,count=1)
    z = np.fromfile(File, dtype=np.float32,count=Np)
    z *= mul[2] 
    if isObservation != 1:
      z += ranges[2][0]
    chk = np.fromfile(File, dtype=np.int32,count=1)
    File.close()

    if isObservation == 1:
      x = x[0:maskIndex]# * 100/300000
      y = y[0:maskIndex]# * 100/300000
      z = z[0:maskIndex]# * 100/300000

    partData = np.column_stack((x,y,z))

    boxLen = mul

    if isObservation == 1:
      props = vp.getSurveyProps(sample.maskFile, sample.zBoundary[0],
                                sample.zBoundary[1], 
                                sample.zBoundary[0], 
                                sample.zBoundary[1], "all",
                                useLCDM=sample.useLCDM)
      boxVol = props[0]
      volNorm = maskIndex/boxVol
    else:
      boxVol = np.prod(boxLen) 
      volNorm = len(x)/boxVol

    isObservationData = isObservation == 1

    return partData, boxLen, volNorm, isObservationData

# -----------------------------------------------------------------------------
def loadPartVel(workDir, sampleDir, sample):
    #print "    Loading particle velocities..."
    sys.stdout.flush()

    infoFile = workDir+"/"+sampleDir+"/zobov_slice_"+sample.fullName+".par"
    File = NetCDFFile(infoFile, 'r')
    isObservation = getattr(File, 'is_observation')

    if isObservation:
      print "No velocities for observations!"
      return -1
   
    vx = File.variables['vel_x'][0:]
    vy = File.variables['vel_y'][0:]
    vz = File.variables['vel_z'][0:]

    File.close()

    partVel = np.column_stack((vx,vy,vz))

    return partVel


# -----------------------------------------------------------------------------
def shiftPart(inPart, newCenter, periodicLine, boxLen):

  part = inPart.copy()

  part[:,0] -= newCenter[0]
  part[:,1] -= newCenter[1]
  part[:,2] -= newCenter[2]

  shiftUs = np.abs(part[:,0]) > boxLen[0]/2.
  if ("x" in periodicLine): part[shiftUs,0] -= \
                            np.copysign(boxLen[0], part[shiftUs,0])
  shiftUs = np.abs(part[:,1]) > boxLen[1]/2.
  if ("y" in periodicLine): part[shiftUs,1] -= \
                            np.copysign(boxLen[1], part[shiftUs,1])
  shiftUs = np.abs(part[:,2]) > boxLen[2]/2.
  if ("z" in periodicLine): part[shiftUs,2] -= \
                            np.copysign(boxLen[2], part[shiftUs,2])

  return part
