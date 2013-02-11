import os

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# CONFIGURATION

# directory for the input simulation/observational particle files
catalogDir = os.getenv("HOME")+"/workspace/Voids/catalogs/multidark/"

# path to HOD code
hodPath = os.getenv("HOME")+"/projects/Voids/hod/HOD.x"

# where to put the final void catalog, figures, and output logs
voidOutputDir = os.getenv("HOME")+"/workspace/Voids/multidark_nocone/"
figDir = os.getenv("PWD")+"/../figs/multidark_nocone/"
logDir = os.getenv("PWD")+"/../logs/multidark_nocone/"

# where to place the pipeline scripts
scriptDir = os.getenv("PWD")+"/multidark_nocone/"

# simulation or observation?
dataType = "simulation"

# available formats for simulation: gadget, multidark
dataFormat = "multidark"
dataUnit = 1 # as multiple of Mpc/h

# place particles on the lightcone?
useLightCone = False

# common filename of particle files
particleFileBase = "mdr1_particles_z"

# list of file numbers for the particle files
# to get particle file name, we take particleFileBase+fileNum
#fileNums = ["0.53"]
fileNums = ["0.0", "0.53", "1.0"]

# redshift of each file in the above list
#redshifts = ["0.53"]
redshifts = ["0.0", "0.53", "1.0"]

# how many independent slices along the z-axis?
numSlices = 1
#numSlices = 4

# how many subdivisions along the x- and y- axis?
#   ( = 2 will make 4 subvolumes for each slice, = 3 will make 9, etc.)
numSubvolumes = 1

# prefix to give all outputs
prefix = "md_"

# list of desired subsamples - these are in unts of h Mpc^-3!
#subSamples = [0.0004]
subSamples = [0.1, 0.05, 0.01, 0.002, 0.001, 0.0004, 0.000175, 0.00001]

# common filename of halo files, leave blank to ignore halos
haloFileBase = "mdr1_halos_z"

# minimum halo mass cuts to apply for the halo catalog
#   use "none" to get all halos
minHaloMasses = [1.2e13]
#minHaloMasses = ["none", 1.2e13]

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
numZobovThreads = 2

# simulation information
numPart = 2048*2048*2048
lbox = 1000 # Mpc/h
omegaM = 0.27
hubble = 0.70

#galDens = 0.000225
hodParmList = [
  {'name'       : "dr9mid", #BOSS: Manera et al. 2012, eq. 26
   'Mmin'       : 0.0,
   'M1'         : 1.e14,
   'sigma_logM' : 0.596,
   'alpha'      : 1.0127,
   'Mcut'       : 1.19399e13,
   'galDens'    : 0.00016,
  },

  {'name'       : "dr7dim2",
   'Mmin'       : 1.99525e12,
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
