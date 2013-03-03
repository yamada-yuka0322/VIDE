#+
#   VIDE -- Void IDEntification pipeline -- ./pipeline/datasets/mock_dr9mid.py
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

# directory for the input simulation/observational particle files
catalogDir = os.getenv("HOME")+"/workspace/Voids/catalogs/mock_dr9mid/"

# path to HOD code
hodPath = os.getenv("HOME")+"/projects/Voids/hod/HOD.x"

# path to mask
maskFile = os.getenv("HOME")+"/workspace/Voids/catalogs/boss/final_boss_mask.fits")

# where to put the final void catalog, figures, and output logs
voidOutputDir = os.getenv("HOME")+"/workspace/Voids/mock_dr9mid/"
figDir = os.getenv("PWD")+"/../figs/mock_dr9mid/"
logDir = os.getenv("PWD")+"/../logs/mock_dr9mid/"

# where to place the pipeline scripts
scriptDir = os.getenv("PWD")+"/mock_dr9mid/"

# simulation or observation?
dataType = "observation"

# available formats for simulation: gadget, multidark
dataFormat = "multidark"
dataUnit = 1 # as multiple of Mpc/h

# place particles on the lightcone?
useLightCone = True

# list of file numbers for the particle files
# to get particle file name, we take particleFileBase+fileNum
fileNums = (("0.53"))

# redshift range of the mock
redshiftRange = (0.53, 0.6)

# prefix to give all outputs
prefix = "mock_"

# common filename of halo files
haloFileBase = "mdr1_halos_z"

# adjust these two parameters given the memory contraints on your system:
#   numZobovDivisions: how many sub-volumes per dimension will zobov process
#   numZobovThreads: how many sub-volumes to process at once?   
numZobovDivisions = 2
numZobovThreads = 2

# simulation information
numPart = 1024*1024*1024
lbox = 1000 # Mpc/h
omegaM = 0.27
hubble = 0.70


# END CONFIGURATION
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
