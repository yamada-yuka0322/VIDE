#!/usr/bin/env python

# prepares input catalogs based on multidark simulations
#   (borrows heavily from generateMock, but doesn't hold much in memory)
# also creates necessary analyzeVoids input files

import numpy as np
import os
import sys
import void_python_tools as vp
import argparse
import imp
import healpy as hp

# ------------------------------------------------------------------------------

def my_import(name):
    mod = __import__(name)
    components = name.split('.')
    for comp in components[1:]:
        mod = getattr(mod, comp)
    return mod

# -----------------------------------------------------------------------------

LIGHT_SPEED = 299792.458

parser = argparse.ArgumentParser(description='options')
parser.add_argument('--scripts', dest='script', action='store_const',
                   const=True, default=False,
                   help='write scripts')
parser.add_argument('--parmFile', dest='parmFile',
                   default="",
                   help='path to parameter file')
args = parser.parse_args()


filename = args.parmFile
print " Loading parameters from", filename
if not os.access(filename, os.F_OK):
  print "  Cannot find parameter file %s!" % filename
  exit(-1)
parms = imp.load_source("name", filename)
globals().update(vars(parms))


#------------------------------------------------------------------------------
def getSampleName(setName, redshift, useVel, iSlice=-1, iVol=-1):

  sampleName = setName

  sampleName += "_z" + redshift

  if iVol != -1: sampleName += "_d" + iVol

  return sampleName

#------------------------------------------------------------------------------
# for given dataset parameters, outputs a script for use with analyzeVoids
def writeScript(setName, dataFileNameBase,
                scriptDir, catalogDir, fileNums, redshifts, numSubvolumes,
                numSlices, useVel, lbox, minRadius, omegaM, subsample=1.0, 
                suffix=".dat"):


  if useVel: setName += "_pv"

  scriptFileName = scriptDir + "/" + setName + ".py"
  scriptFile = open(scriptFileName, 'w')

  scriptFile.write("""#!/usr/bin/env/python
import os
from void_python_tools.backend.classes import *

continueRun = True # set to True to enable restarting aborted jobs
startCatalogStage = 1
endCatalogStage   = 4
               
startAPStage = 1
endAPStage = 7

ZOBOV_PATH = os.getenv("PWD")+"/../zobov/"
CTOOLS_PATH = os.getenv("PWD")+"/../c_tools/"
freshStack = True
errorBars = "CALCULATED"
numIncoherentRuns = 100
ranSeed = 101010
useLCDM = False
bias = 1.16

dataPortions = ["central"]
dataSampleList = []
""")


  dataInfo = """
setName = "{setName}"

workDir = "{voidOutputDir}/{setName}/"
inputDataDir = "{inputDataDir}"
figDir = "{figDir}/{setName}/"
logDir = "{logDir}/{setName}/"

numZobovDivisions = {numZobovDivisions}
numZobovThreads = {numZobovThreads}
               """
  scriptFile.write(dataInfo.format(setName=setName, figDir=figDir,
                                   logDir=logDir, voidOutputDir=voidOutputDir,
                                   inputDataDir=catalogDir, 
                                   numZobovDivisions=numZobovDivisions,
                                   numZobovThreads=numZobovThreads))

  sampleInfo = """
newSample = Sample(dataFile = "{dataFile}",
                   dataFormat = "{dataFormat}",
                   dataUnit = {dataUnit},
                   fullName = "{sampleName}",
                   nickName = "{sampleName}",
                   dataType = "simulation",
                   zBoundary = ({zMin}, {zMax}),
                   zRange    = ({zMin}, {zMax}),
                   zBoundaryMpc = ({zMinMpc}, {zMaxMpc}),
                   omegaM    = {omegaM},
                   minVoidRadius = {minRadius},
                   includeInHubble = True,
                   partOfCombo = False,
                   isCombo = False,
                   boxLen = {boxLen},
                   usePecVel = {usePecVel},
                   numSubvolumes = {numSubvolumes},
                   mySubvolume = "{mySubvolume}",
                   useLightCone = {useLightCone},
                   subsample = {subsample})
dataSampleList.append(newSample)
newSample.addStack({zMin}, {zMax}, {minRadius}  , {minRadius}+2, True, False)
newSample.addStack({zMin}, {zMax}, {minRadius}  , {minRadius}+4, True, False)
newSample.addStack({zMin}, {zMax}, {minRadius}+2, {minRadius}+6, True, False)
newSample.addStack({zMin}, {zMax}, {minRadius}+6, {minRadius}+10, True, False)
newSample.addStack({zMin}, {zMax}, {minRadius}+10, {minRadius}+18, True, False)
newSample.addStack({zMin}, {zMax}, {minRadius}+18, {minRadius}+24, True, False)
               """
  for (iFile, redshift) in enumerate(redshifts):
    fileNum = fileNums[iFile]
    zBox = float(redshift)
    Om = float(omegaM)
    zBoxMpc = LIGHT_SPEED/100.*vp.angularDiameter(zBox, Om=Om)
    boxMaxMpc = zBoxMpc + lbox

    # converter from redshift to comoving distance   
    zVsDY = np.linspace(0., zBox+8*lbox*100./LIGHT_SPEED, 10000)
    zVsDX = np.zeros(len(zVsDY))
    for i in xrange(len(zVsDY)):
      zVsDX[i] = vp.angularDiameter(zVsDY[i], Om=Om)

    if useLightCone:
      boxWidthZ = np.interp(vp.angularDiameter(zBox,Om=Om)+100. / \
                  LIGHT_SPEED*lbox, zVsDX, zVsDY)-zBox
      dzSafe = 0.03
    else:
      boxWidthZ = lbox*100./LIGHT_SPEED
      dzSafe = 0.0

    for iSlice in xrange(numSlices):
      sliceMin = zBox + dzSafe + iSlice*(boxWidthZ-dzSafe)/numSlices
      sliceMax = zBox + dzSafe + (iSlice+1)*(boxWidthZ-dzSafe)/numSlices

      sliceMinMpc = sliceMin*LIGHT_SPEED/100.
      sliceMaxMpc = sliceMax*LIGHT_SPEED/100.

      sliceMin = "%0.2f" % sliceMin
      sliceMax = "%0.2f" % sliceMax
      sliceMinMpc = "%0.1f" % sliceMinMpc
      sliceMaxMpc = "%0.1f" % sliceMaxMpc

      dataFileName = dataFileNameBase + fileNum + suffix

      for iX in xrange(numSubvolumes):
        for iY in xrange(numSubvolumes):

          mySubvolume = "%d%d" % (iX, iY)

          sampleName = getSampleName(setName, sliceMin, useVel,
                                     iSlice=iSlice, iVol=mySubvolume)

          scriptFile.write(sampleInfo.format(dataFile=dataFileName,
                                         dataFormat=dataFormat,
                                         dataUnit=dataUnit,
                                         sampleName=sampleName,
                                         zMin=sliceMin,
                                         zMax=sliceMax,
                                         zMinMpc=sliceMinMpc,
                                         zMaxMpc=sliceMaxMpc,
                                         omegaM=Om,
                                         boxLen=lbox,
                                         usePecVel=useVel,
                                         minRadius=minRadius,
                                         numSubvolumes=numSubvolumes,
                                         mySubvolume=mySubvolume,
                                         useLightCone=useLightCone,
                                         subsample=subsample))

  scriptFile.close()
  return


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
if not os.access(scriptDir, os.F_OK): os.mkdir(scriptDir)
if not os.access(catalogDir, os.F_OK): os.mkdir(catalogDir)

# -----------------------------------------------------------------------------
# now the SDSS HOD
parFileText = """
% cosmology
OMEGA_M {omegaM}
HUBBLE  {hubble}
OMEGA_B 0.0469
SIGMA_8 0.82
SPECTRAL_INDX 0.95
ITRANS 5
REDSHIFT {redshift}

% halo definition
%DELTA_HALO 200
DELTA_HALO 740.74 % 200/Om_m
M_max 1.00E+16

% fit function types
pdfs 11
pdfc 2
EXCLUSION 4

% hod parameters 
M_min {Mmin}
GALAXY_DENSITY 0.0111134 % computed automatically if M_min set, use for sanity
M1 {M1}   
sigma_logM {sigma_logM}
alpha {alpha}
M_cut {Mcut}

% simulation info
real_space_xi 1
HOD 1
populate_sim 1
HaloFile {haloFile}
RESOLUTION {numPartPerSide}
BOX_SIZE   {boxSize}

% output
root_filename hod
               """

print " Doing DR9 HOD"

    # these parameters come from Manera et al 2012, eq. 26
parFileName = "./hod.par"
parFile = open(parFileName, 'w')
haloFile = catalogDir+haloFileBase+fileNums[iRedshift]
parFile.write(parFileText.format(omegaM=omegaM,
                                     hubble=hubble,
                                     redshift=redshift,
                                     Mmin=1.23e13,
                                     M1=1.e14,
                                     sigma_logM=0.596,
                                     alpha=1.0127,
                                     Mcut=1.19399e13,
                                     haloFile=haloFile,
                                     numPartPerSide=numPart**(1/3.),
                                     boxSize=lbox))
parFile.close()

os.system(hodPath+" "+parFileName+">& /dev/null")

# now place these particles on a lightcone, restrict redshift range, apply mask
mask = hp.read_map(maskFile)
nside = hp.get_nside(mask)

inFile  = open('hod.mock', 'r') 
outFile = open(catalogDir+"/mock.out"))

zBox = float(redshiftRange[0])
Om = float(omegaM)

# converter from redshift to comoving distance   
zVsDY = np.linspace(0., zBox+8*lbox*100./LIGHT_SPEED, 10000)
zVsDX = np.zeros(len(zVsDY))
for i in xrange(len(zVsDY)):
  zVsDX[i] = vp.angularDiameter(zVsDY[i], Om=Om)

for line in inFile:
  line = line.split(',')
  x  = float(line[0]) - lbox/2.
  y  = float(line[1]) - lbox/2.
  z  = float(line[2]) - lbox/2.
  vz = float(line[5])

  zBoxInMpc = vp.angularDiameter(zBox, Om=Om)

  redshift = np.sqrt(x*x + y*y + z*z) 
  redshift = np.interp(zBoxInMpc+100./LIGHT_SPEED*redshift, zVsDX, zVsDY)

  if redshift < redshiftRange[0] or redshift > redshiftRange[1]: continue

  RA = np.atan((y-lbox/2.)/(x-lbox/2.)) * 100/np.pi + 180.
  Dec = np.asin((z-lboc/2.)/(redshift*LIGHT_SPEED/100.)) * 180/np.pi

  phi = np.pi/180. * RA
  theta = np.pi/2. - Dec*np.pi/180.
  pos = np.zeros((3))

  pix = hp.ang2pix(nside, theta, phi)
  if mask[pix] <= 0: continue   
  
  print >> outFile, RA, Dec, redshift*LIGHT_SPEED, 0., x, y, z
 
inFile.close()
outFile.close()

os.system("rm ./hod.*")


