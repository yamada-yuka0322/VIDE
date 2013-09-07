#!/usr/bin/env python
#+
#   VIDE -- Void IDEntification pipeline -- ./pipeline/apAnalysis.py
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

# plots radial density profiles centered on baseVoid.
#   requires makeCocenteredProfiles to be run first!

import imp
import pickle
import os
import numpy as np
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from void_python_tools.backend import *
from util import *
from globalOptions import *
from scipy.optimize import curve_fit

matplotlib.rcParams.update({'font.size': 20})

# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='Analyze.')
parser.add_argument('--parm', dest='parm', default='datasetsToAnalyze.py', help='path to parameter file')
parser.add_argument('--show', dest='showPlot', action='store_const',
                   const=True, default=False,
                   help='display the plot (default: just write eps)')
args = parser.parse_args()

# -----------------------------------------------------------------------------
# Lavaux & Wandelt (2012) profile
def LWProfile(r, A0, A3, alpha):
  return A0 + A3*r**3.
  #return A0 + A3*r**alpha

# -----------------------------------------------------------------------------
# http://arxiv.org/pdf/astro-ph/0508297v1.pdf eq. 5
def PadillaProfile(r, A1, A2, alpha):
  return 1.5-A1*np.exp(-(A2*r)**alpha)


# -----------------------------------------------------------------------------
# plot a slice of the density around the void in baseIDList, 
#   with any voids in the slice shown and any voids in baseIDList flagged
def plotProfiles(baseSample, stack, sampleList,
                 figDir, showPlot, outputDir, doTheory):

  thisRadius = str(stack.rMin) + "-" + str(stack.rMax)
  plotName = "1dprofile_cocenter_" + baseSample.fullName+"_"+thisRadius

  filename = "1dprofile_cocenter_" + baseSample.fullName+"_"+thisRadius
  npzfile = np.load(outputDir+"/1dprofile_cocentered_"+plotName+".dat.npz")
  profileList = npzfile['arr_0']
  radii = npzfile['arr_1']

  plt.clf()

  plt.xlabel(r"$R/R_{eff}$")
  #plt.xlabel(r"$R/R_{v,\mathrm{max}}$")
  plt.ylabel(r"$n / \bar n$")
  plt.xlim(xmin=0.0, xmax=2.5)
  plt.ylim(ymin=0.0, ymax=1.7)
  #plt.xscale('log')

  for (iSample, sample) in enumerate(sampleList):
    lineTitle = sample.nickName[:-10]
    if "DM LowDen" in lineTitle or "DM HighDen" in lineTitle: continue

    thisPlotTitle = r"$R_{eff}$ = "+thisRadius+ r" $h^{-1}$Mpc"
    legendTitle = "Fixed Center Samples"

    thisProfile = profileList[iSample]
    if np.all(thisProfile == 0.): continue

    rV = (stack.rMin + stack.rMax)/2.

    if len(radii) > 0:
      scaledRadii = radii/rV
      plt.plot(scaledRadii, thisProfile, label=lineTitle, 
               color=colorList[iSample],
               linewidth=linewidth)

      if doTheory:
        nBins = len(scaledRadii)/2
        try:
          popt, pcov = curve_fit(PadillaProfile, scaledRadii[0:nBins],
                             thisProfile[0:nBins], maxfev=10000, xtol=5.e-3)
        except RuntimeError:
          print "Warning: no convergence reached"

        label = r"A1=%.2f, A2=%.2f, $\alpha$=%.2f" % (popt[0], popt[1], popt[2])
        rho = PadillaProfile(scaledRadii, popt[0], popt[1], popt[2])
        plt.plot(scaledRadii, rho, '--', label=label,
                 color=colorList[iSample], linewidth=2)


  plt.title(thisPlotTitle + " (center from " + plotTitle + ")", fontsize=18)
  plt.legend(title = legendTitle, loc = "lower right", prop={'size':16})

  plt.savefig(figDir+"/fig_"+plotName+".pdf", bbox_inches="tight")
  plt.savefig(figDir+"/fig_"+plotName+".eps", bbox_inches="tight")
  plt.savefig(figDir+"/fig_"+plotName+".png", bbox_inches="tight")

  if (showPlot):
    os.system("display %s" % figDir+"/fig_"+plotName+".png")

  return

# -----------------------------------------------------------------------------

filename = args.parm
print " Loading parameters from", filename
if not os.access(filename, os.F_OK):
  print "  Cannot find parameter file %s!" % filename
  exit(-1)
parms = imp.load_source("name", filename)
globals().update(vars(parms))

if not os.access(outputDir, os.F_OK):
  os.makedirs(outputDir)

if not os.access(logDir, os.F_OK):
  os.makedirs(logDir)

if not os.access(figDir, os.F_OK):
  os.makedirs(figDir)

# get list of base voids
with open(workDir+baseSampleDir+"/sample_info.dat", 'rb') as input:
  baseSample = pickle.load(input)
  baseSampleName = baseSample.fullName
  baseVoidList = np.loadtxt(workDir+baseSampleDir+"/centers_central_"+\
                    baseSampleName+".out")

sampleList = []
for sampleDir in sampleDirList:
  if compareSampleTag in sampleDir: continue
  with open(workDir+sampleDir+"/sample_info.dat", 'rb') as input:
    sampleList.append(pickle.load(input))

sampleDirList.insert(0,baseSampleDir)
sampleList.insert(0,baseSample)

# pick our void sample
for stack in baseSample.stacks:
  print " Stack:", stack.rMin, "-", stack.rMax

  # plot these profiles
  print "  Plotting..."
  sys.stdout.flush()
  plotProfiles(baseSample, stack, sampleList,
               figDir, args.showPlot, outputDir, doTheory)
    
print " Done!"

