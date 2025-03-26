#!/usr/bin/env python
#+
#   VIDE -- Void IDentification and Examination -- ./python_tools/vide/voidUtil/matchUtil.py
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
__all__=['compareCatalogs',]

from vide.backend import *
import imp
import pickle
import os
import matplotlib.pyplot as plt
import numpy as np
import argparse

def compareCatalogs(baseCatalogDir, compareCatalogDir, 
                    outputDir="./", logDir="./", 
                    matchMethod="useID", dataPortion="central",
                    strictMatch=True, 
                    overlapFrac=0.25,
                    pathToCTools="../../../c_tools"):

# reports the overlap between two void catalogs
#  baseCatalogDir: directory of catalog 1
#  compareCatagDir: directory of catalog 2
#  outputDir: directory to place outputs
#  logDir: directory to place log files
#  matchMethod: "useID" to use unique IDs, "prox" to use overlap of Voronoi cells
#  dataPortion: "central" or "all"
#  strictMatch: if True, only attempt to match to trimmed catalog
#  overlapFrac: threshold fraction of Voronoi radius to count as matched
#  pathToCTools: path to location of VIDE c_tools directory

  if not os.access(outputDir, os.F_OK):
    os.makedirs(outputDir)
  
  if not os.access(logDir, os.F_OK):
    os.makedirs(logDir)
  
  outFileName = outputDir + "/" + "voidOverlap" #+ ".dat"
  
  with open(baseCatalogDir+"/sample_info.dat", 'rb') as input:
    baseSample = pickle.load(input)
  
  with open(compareCatalogDir+"/sample_info.dat", 'rb') as input:
    sample = pickle.load(input)
  
  print(" Comparing", baseSample.fullName, "to", sample.fullName, "...", end=' ')
  sys.stdout.flush()
  
  sampleName = sample.fullName
  
  binPath = pathToCTools+"/analysis/voidOverlap"
  logFile = logDir+"/compare_"+baseSample.fullName+"_"+sampleName+".out"
  stepOutputFileName = outFileName + "_" + baseSample.fullName + "_" + \
                       sampleName + "_"
  
  launchVoidOverlap(baseSample, sample, baseCatalogDir, 
                    compareCatalogDir, binPath, 
                    thisDataPortion=dataPortion, logFile=logFile,
                    continueRun=False, 
                    outputFile=stepOutputFileName,
                    matchMethod=matchMethod,
                    overlapFrac=overlapFrac,
                    strictMatch=strictMatch)
  
  print(" Done!")
  return
