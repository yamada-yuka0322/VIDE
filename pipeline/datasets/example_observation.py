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

# a global name to give 
#catalogName = "lcdm"


# directory for input data files
inputDataDir = os.getenv("HOME")+"/workspace/Voids/catalogs/nyuvagc/"

# void catalog output directory
workDir      = os.getenv("HOME")+"/workspace/Voids/sdss_dr7LCDM/"

# output directory for log files
logDir = os.getenv("PWD")+"/../logs/sdss_dr7LCDM"

# output directory for figures
figDir = os.getenv("PWD")+"/../figs/sdss_dr7LCDM"

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
                   dataFile = "filename.dat"

                   # full name for this sample
                   fullName = "lss.dr72dim1.dat",

                   # a convenient nickname
                   nickName = "dim1",

                   # don't change this
                   dataType = "observation",

                   # HEALpix mask file
                   maskFile = inputDataDir+"/healpix/rast_window_512.fits",

                   # max and min redshifts of galaxies in your sample
                   zBoundary = (0.0, 0.05),

                   # max and min redshifts where you want to find voids
                   zRange = (0.0, 0.05),
 
                   # TODO 
                   skyFraction = 0.19,

                   # leave this at -1 for mean particle separation, or 
                   # specify your own in Mpc/h
                   minVoidRadius = -1,

                   # density of mock particles in cubic Mpc/h
                   fakeDensity = 0.01,

                   )
dataSampleList.append(newSample)

# repeat the above block for any other samples
