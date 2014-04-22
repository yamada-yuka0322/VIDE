#+
#   VIDE -- Void IDEntification pipeline -- ./pipeline/datasets/mergertree.py
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
import os

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# CONFIGURATION

startCatalogStage = 0
endCatalogStage = 0
continueRun = False

startAPStage = 1
endAPStage = 2

# directory for the input simulation/observational particle files
catalogDir = os.getenv("HOME")+"/workspace/Voids/catalogs/mergertree1024/"

# path to HOD code
hodPath = os.getenv("HOME")+"/projects/Voids/hod/HOD.x"

# where to put the final void catalog, figures, and output logs
voidOutputDir = os.getenv("HOME")+"/workspace/Voids/mergertree1024/"
figDir = os.getenv("PWD")+"/../figs/mergertree1024/"
logDir = os.getenv("PWD")+"/../logs/mergertree1024/"

# where to place the pipeline scripts
scriptDir = os.getenv("PWD")+"/mergertree1024/"

# simulation or observation?
dataType = "simulation"

# available formats for simulation: gadget, mergertree
dataFormat = "sdf"
dataUnit = 1 # as multiple of Mpc/h

# place particles on the lightcone?
useLightCone = False

# also do peculiar velocities?
doPecVel = False

# common filename of particle files
particleFileBase = "mf_4s_1G_1k_NNNNN"
particleFileDummy = 'NNNNN'

# list of file numbers for the particle files
# to get particle file name, we take particleFileBase+fileNum
fileNums = ["1.000"]

# redshift of each file in the above list
redshifts = ["0.0"]

# how many independent slices along the z-axis?
numSlices = 1
#numSlices = 4
numAPSlices = 1

# how many subdivisions along the x- and y- axis?
#   ( = 2 will make 4 subvolumes for each slice, = 3 will make 9, etc.)
numSubvolumes = 1

# prefix to give all outputs
prefix = "mt_"

# list of desired subsamples - these are in unts of h Mpc^-3!
subSamples = [0.1, 0.05, 0.02, 0.01, 0.004, 0.002, 0.001, 0.0003, 0.0001]
#doSubSampling = False # do the subsampling in preparation script?
doSubSampling = True # do the subsampling in preparation script?

# common filename of halo files, leave blank to ignore halos
haloFileBase = "mf_4s_1G_1k_bgc2_NNNNN.sdf"
haloFileDummy = 'NNNNN'

# minimum halo mass cuts to apply for the halo catalog
#   use "none" to get all halos
minHaloMasses = ["none", 1.2e13]
#minHaloMasses = [7.0e11, 1.2e13]

# locations of data in the halo catalog
haloFileMCol  = 6
haloFileXCol  = 0
haloFileYCol  = 1
haloFileZCol  = 2
haloFileVXCol = 3
haloFileVYCol = 4
haloFileVZCol = 5
haloFileColSep = ','
haloFileNumComLines = 0

# adjust these two parameters given the memory contraints on your system:
#   numZobovDivisions: how many sub-volumes per dimension will zobov process
#   numZobovThreads: how many sub-volumes to process at once?   
numZobovDivisions = 4
numZobovThreads = 4

# simulation information
numPart = 1024*1024*1024
lbox = 999.983 # Mpc/h
omegaM = 0.2847979853038958
hubble = 0.6962

#galDens = 0.000225
hodParmList = [
  {'name'       : "LowRes", #BOSS: Manera et al. 2012, eq. 26
   'Mmin'       : 0.0,
   'M1'         : 1.e14,
   'sigma_logM' : 0.596,
   'alpha'      : 1.0127,
   'Mcut'       : 1.19399e13,
   'galDens'    : 0.0002,
  },

  {'name'       : "HighRes",
   'Mmin'       : 0.0,
   #'Mmin'       : 1.99525e12,
   'M1'         : 3.80189e13,
   'sigma_logM' : 0.21,
   'alpha'      : 1.12,
   'Mcut'       : 6.91831e11,
   'galDens'    : 0.02,
  }
]

# END CONFIGURATION
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
