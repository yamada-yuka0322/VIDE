#+
#   VIDE -- Void IDentification and Examination -- ./pipeline/datasets/example_observation.py
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

import os
import numpy as np
from void_python_tools.backend.classes import *

# if True, will scan log files for last known completed state and run from there
continueRun = False

# stages:
#   1 : extract redshift slices from data
#   2 : void extraction using zobov
#   3 : removal of small voids and voids near the edge 
startCatalogStage = 1
endCatalogStage   = 3

# directory for input data files
inputDataDir = os.getenv("PWD")+"/../examples/"

# void catalog output directory
workDir = os.getenv("PWD")+"/../examples/example_observation/"

# output directory for log files
logDir = os.getenv("PWD")+"/../logs/example_observation/"

# output directory for figures
figDir = os.getenv("PWD")+"/../figs/example_observation/"

# you need to set these manually: point to ZOBOV and C_TOOLS in VIDE directory
ZOBOV_PATH = os.getenv("PWD")+"/../zobov/"
CTOOLS_PATH = os.getenv("PWD")+"/../c_tools/"

# if true, convert to comoving space using LCDM cosmology
useComoving = True

# optimization: maximum number of parallel threads to use
numZobovThreads = 2

# optimization: number of subdivisions of the box
numZobovDivisions = 2

# don't change this
dataSampleList = []

# define your volume-limited samples
newSample = Sample(
                   # path to galaxy file is inputDataDir+dataFile
                   dataFile = "example_observation.dat",

                   # full name for this sample
                   fullName = "example_observation",

                   # a convenient nickname
                   nickName = "exobs",

                   # don't change this
                   dataType = "observation",

                   # assume sample is volume-limited?
                   volumeLimited = True,

                   # HEALpix mask file
                   maskFile = inputDataDir+"/example_observation_mask.fits",

                   # radial selection function (if not volume limited)
                   selFunFile = None,

                   # max and min redshifts of galaxies in your sample
                   zBoundary = (0.0, 0.15),

                   # max and min redshifts where you want to find voids
                   zRange = (0.1, 0.15),
 
                   # leave this at -1 for mean particle separation, 
                   # or specify your own in Mpc/h
                   minVoidRadius = -1,

                   # density of mock particles in cubic Mpc/h
                   # (make this as high as you can afford)
                   fakeDensity = 0.05,

                   )
dataSampleList.append(newSample)

# repeat the above block for any other samples
