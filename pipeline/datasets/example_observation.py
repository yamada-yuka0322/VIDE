#+
#   VIDE -- Void IDentification and Examination -- ./pipeline/sdss/sdss_dr7LCDM.py
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
#!/usr/bin/env python

# does analysis in real space assuming LCDM cosmology

import os
import numpy as np
from void_python_tools.backend.classes import *

# does full void analysis. Also generates 2d/1d stacked plots and hubble diagram

# to combine multiple samples:
#   For example, combining dim from 0.0-0.1 and bright from 0.1-0.2:
#     1) add dim+bright to sampleList - "+" is required to detect combo!

# if True, will scan log files for last known completed state and run from there
continueRun = False

# stages:
#   1 : extract redshift slices from data
#   2 : void extraction using zobov
#   3 : removal of small voids and voids near the edge 
#   4 : void stacking and combining
#   5 : 1-d and 2-d profile generation
#   6 : fitting and plotting stacked voids
#   7 : determination of error bars (if requested)
#   8 : hubble diagram generation
#   9 : liklihood determination
startCatalogStage = 1
endCatalogStage   = 1

catalogName = "lcdm"

# where main data files are present/will go
workDir      = os.getenv("HOME")+"/workspace/Voids/sdss_dr7LCDM/"
inputDataDir = os.getenv("HOME")+"/workspace/Voids/catalogs/nyuvagc/"

logDir = os.getenv("PWD")+"/../logs/sdss_dr7LCDM"
figDir = os.getenv("PWD")+"/../figs/sdss_dr7LCDM"

ZOBOV_PATH = os.getenv("PWD")+"/../zobov/"
CTOOLS_PATH = os.getenv("PWD")+"/../c_tools/"

freshStack = True # "no" to continue building onto existing stacked voids

errorBars = "CALCULATED" # "ESTIMATED" for error bars from multiple null tests
#errorBars = "ESTIMATED" # "ESTIMATED" for error bars from multiple null tests
                        # "CALCULATED" to use fitting procedure 
numIncoherentRuns = 100

ranSeed = 101010

useComoving = True # if true, convert to real space using LCDM cosmology

# try different portions of data. available options: "all, central, edge"
dataPortions = ["central", "all"]

numZobovDivisions = 2
numZobovThreads = 2

dataSampleList = []


newSample = Sample(fullName = "lss.dr72dim1.dat",
                   nickName = "dim1",
                   dataType = "observation",
                   maskFile = inputDataDir+"/healpix/rast_window_512.fits",
                   selFunFile = inputDataDir+"/czselfunc.all.dr72dim1.dat",
                   zBoundary = (0.0, 0.05),
                   zRange = (0.0, 0.05),
                   skyFraction = 0.19,
                   minVoidRadius = 3.479,
                   fakeDensity = 0.01,
                   volumeLimited = True,
                   includeInHubble = False,
                   partOfCombo = True,
                   isCombo = False,
                   comboList= None)
newSample.addStack(0.0, 0.1, 8, 12, False, True)
newSample.addStack(0.0, 0.1, 12, 16, False, True)
newSample.addStack(0.0, 0.1, 16, 20, False, True)
newSample.addStack(0.0, 0.1, 20, 28, False, True)
dataSampleList.append(newSample)

newSample = Sample(fullName = "lss.dr72dim2.dat",
                   nickName = "dim2",
                   dataType = "observation",
                   maskFile = inputDataDir+"/healpix/rast_window_512.fits",
                   selFunFile = inputDataDir+"/czselfunc.all.dr72dim2.dat",
                   zBoundary = (0.0, 0.1),
                   zRange = (0.05, 0.1),
                   skyFraction = 0.19,
                   minVoidRadius = 4.723,
                   fakeDensity = 0.01,
                   volumeLimited = True,
                   includeInHubble = False,
                   partOfCombo = True,
                   isCombo = False,
                   comboList= None)
newSample.addStack(0.0, 0.1, 8, 12, False, True)
newSample.addStack(0.0, 0.1, 12, 16, False, True)
newSample.addStack(0.0, 0.1, 16, 20, False, True)
newSample.addStack(0.0, 0.1, 20, 28, False, True)
dataSampleList.append(newSample)

newSample = Sample(fullName = "lss.dr72dim1+dim2.dat",
                   nickName = "dim1+dim2",
                   dataType = "observation",
                   maskFile = inputDataDir+"/healpix/rast_window_512.fits",
                   selFunFile = inputDataDir+"/czselfunc.all.dr72dim.dat",
                   zBoundary = (0.0, 0.1),
                   zRange = (0.0, 0.1),
                   skyFraction = 0.19,
                   minVoidRadius = 4.723,
                   fakeDensity = 0.01,
                   profileBinSize = "auto",
                   volumeLimited = True,
                   includeInHubble = True,
                   partOfCombo = False,
                   isCombo = True,
                   comboList= ("lss.dr72dim1.dat", "lss.dr72dim2.dat"))
newSample.addStack(0.0, 0.1, 8, 12, True, False)
newSample.addStack(0.0, 0.1, 12, 16, True, False)
newSample.addStack(0.0, 0.1, 16, 20, True, False)
newSample.addStack(0.0, 0.1, 20, 28, True, False)
dataSampleList.append(newSample)


newSample = Sample(fullName = "lss.dr72bright1.dat",
                   nickName = "bright1",
                   dataType = "observation",
                   maskFile = inputDataDir+"/healpix/rast_window_512.fits",
                   selFunFile = inputDataDir+"/czselfunc.all.dr72bright1.dat",
                   zBoundary = (0.0, 0.15),
                   zRange = (0.1, 0.15),
                   skyFraction = 0.19,
                   minVoidRadius = 6.763,
                   fakeDensity = 0.01,
                   volumeLimited = True,
                   includeInHubble = False,
                   partOfCombo = True,
                   isCombo = False,
                   comboList= None)
newSample.addStack(0.1, 0.2, 13, 20, False, True)
newSample.addStack(0.1, 0.2, 20, 28, False, True)
newSample.addStack(0.1, 0.2, 28, 36, False, True)
newSample.addStack(0.1, 0.2, 36, 44, False, True)
dataSampleList.append(newSample)

newSample = Sample(fullName = "lss.dr72bright2.dat",
                   nickName = "bright2",
                   dataType = "observation",
                   maskFile = inputDataDir+"/healpix/rast_window_512.fits",
                   selFunFile = inputDataDir+"/czselfunc.all.dr72bright2.dat",
                   zBoundary = (0.0, 0.2),
                   zRange = (0.15, 0.2),
                   skyFraction = 0.19,
                   minVoidRadius = 10.844,
                   fakeDensity = 0.001,
                   volumeLimited = True,
                   includeInHubble = False,
                   partOfCombo = True,
                   isCombo = False,
                   comboList= None)
newSample.addStack(0.1, 0.2, 13, 20, False, True)
newSample.addStack(0.1, 0.2, 20, 28, False, True)
newSample.addStack(0.1, 0.2, 28, 36, False, True)
newSample.addStack(0.1, 0.2, 36, 44, False, True)
dataSampleList.append(newSample)

newSample = Sample(fullName = "lss.dr72bright1+bright2.dat",
                   nickName = "bright1+bright2",
                   dataType = "observation",
                   maskFile = inputDataDir+"/healpix/rast_window_512.fits",
                   zBoundary = (0.0, 0.2),
                   zRange = (0.1, 0.2),
                   skyFraction = 0.19,
                   minVoidRadius = -1,
                   fakeDensity = -1,
                   profileBinSize = "auto",
                   volumeLimited = True,
                   includeInHubble = True,
                   partOfCombo = False,
                   isCombo = True,
                   comboList= ("lss.dr72bright1.dat", "lss.dr72bright2.dat"))
newSample.addStack(0.1, 0.2, 13, 20, True, False)
newSample.addStack(0.1, 0.2, 20, 28, True, False)
newSample.addStack(0.1, 0.2, 28, 36, True, False)
newSample.addStack(0.1, 0.2, 36, 44, True, False)
dataSampleList.append(newSample)

newSample = Sample(fullName = "lss.dr72lrgdim.dat",
                   nickName = "lrgdim",
                   dataType = "observation",
                   maskFile = inputDataDir+"/healpix/rast_window_512.fits",
                   selFunFile = inputDataDir+"/czselfunc.all.dr72lrgdim.dat",
                   zBoundary = (0.16, 0.36),
                   zRange = (0.16, 0.36),
                   skyFraction = 0.19,
                   minVoidRadius = 24,
                   fakeDensity = 0.001, 
                   profileBinSize = 2, # KEEP
                   volumeLimited = True,
                   includeInHubble = True,
                   partOfCombo = False,
                   isCombo = False,
                   comboList= None)
newSample.addStack(0.16, 0.36, 54, 70, True, False, zMinPlot=0.2)
newSample.addStack(0.16, 0.36, 70, 106, True, False, zMinPlot=0.2)
dataSampleList.append(newSample)

newSample = Sample(fullName = "lss.dr72lrgbright.dat",
                   nickName = "lrgbright",
                   dataType = "observation",
                   maskFile = inputDataDir+"/healpix/rast_window_512.fits",
                   selFunFile = inputDataDir+"/czselfunc.all.dr72lrgbright.dat",
                   zBoundary = (0.36, 0.44),
                   zRange = (0.36, 0.44),
                   skyFraction = 0.19,
                   minVoidRadius = 38,
                   fakeDensity = 0.001,
                   profileBinSize = "auto",
                   volumeLimited = True,
                   includeInHubble = False,
                   partOfCombo = False,
                   isCombo = False,
                   comboList= None)
newSample.addStack(0.36, 0.44, 76, 92, False, False)
newSample.addStack(0.36, 0.44, 92, 108, False, False)
dataSampleList.append(newSample)


