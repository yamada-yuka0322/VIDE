#+
#   VIDE -- Void IDentification and Examination -- ./pipeline/datasets/mergertree.py
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
# DEFAULT CONFIGURATION

datasetName = ""

startCatalogStage = 1
endCatalogStage = 3

startAPStage = 1
endAPStage = 5

continueRun = True
dataPortions = ["central"]

# auto for AP analysis, "fixed" for catalog comparison
stackMode = "fixed"

# directory for the input simulation/observational particle files
catalogDir = os.getenv("HOME")+"/workspace/Voids/catalog/"

# path to HOD code
hodPath = os.getenv("HOME")+"/projects/Voids/hod/HOD.x"

# where to put the final void catalog, figures, and output logs
voidOutputDir = os.getenv("HOME")+"/workspace/Voids//"
figDir = os.getenv("PWD")+"/../figs/"
logDir = os.getenv("PWD")+"/../logs/"

# where to place the pipeline scripts
scriptDir = os.getenv("PWD")+"/scripts//"

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
# use a placeholder (such as NNNNN as shown below) to select the different
#   filenames. For example, if we have partFile01, partFile02, etc., 
#   then particleFileBase = 'partFileNN'
#        particleFileDummy = 'NN'
#        fileNums = ["01", "02"] 
particleFileBase = "mf_4s_1G_512_NNNNN"
particleFileDummy = 'NNNNN'

# list of file numbers for the particle files
# to get particle file name, we replace particleFileDummy with fileNum
fileNums = ["0.667", "0.500"]

# redshift of each file in the above list
redshifts = ["0.5", "1.0"]

# how many independent slices along the z-axis?
numSlices = 1

# how many slices for analysis?
numAPSlices = 1

# how many subdivisions along the x- and y- axis?
#   ( = 2 will make 4 subvolumes for each slice, = 3 will make 9, etc.)
numSubvolumes = 1

# prefix to give all outputs
prefix = "mt_"

# list of desired subsamples - see subSamplingMode parameter
subSamples = [1.0]

doSubSampling = True # do the subsampling in preparation script?
                     # if False, generateMock will do the subsampling

# if 'absolute', subSamples are given in particles per cubic Mpc/h
# if 'relative', subSamples are given as a fraction of input particles
subSampleMode = "absolute" 

# common filename of halo files, leave blank to ignore halos
haloFileBase = "mf_4s_1G_512_bgc2_NNNNN.sdf"
haloFileDummy = 'NNNNN'

# minimum halo mass cuts to apply for the halo catalog
#   use "none" to get all halos
minHaloMasses = [1.2e13]

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
haloFilePosRescale = 1.0 # rescaling necessary to get Mpc/h

# adjust these two parameters given the memory contraints on your system:
#   numZobovDivisions: how many sub-volumes per dimension will zobov process
#   numZobovThreads: how many sub-volumes to process at once?   
numZobovDivisions = 2
numZobovThreads = 2

# simulation information
numPart = 512*512*512
lbox = 999.983 # Mpc/h
omegaM = 0.2847979853038958
hubble = 0.6962

#galDens = 0.000225
hodParmList = [
  {'name'       : "dr9mid", #BOSS: Manera et al. 2012, eq. 26
   'Mmin'       : 0.0,
   'M1'         : 1.e14,
   'sigma_logM' : 0.596,
   'alpha'      : 1.0127,
   'Mcut'       : 1.19399e13,
   'galDens'    : 0.0002,
  },
]

# END CONFIGURATION
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
