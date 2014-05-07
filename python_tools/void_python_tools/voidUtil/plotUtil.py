#+
#   VIDE -- Void IDentification and Examination -- ./python_tools/void_python_tools/plotting/plotTools.py
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
__all__=['plotRedshiftDistribution', 'plotSizeDistribution', 'plot1dProfiles',
         'plotMarg1d', 'plotNumberDistribution', 'plotVoidDistribution']

from void_python_tools.backend.classes import *
from plotDefs import *
import numpy as np
import os
import pylab as plt
import void_python_tools.apTools as vp

def fill_between(x, y1, y2=0, ax=None, **kwargs):
    """Plot filled region between `y1` and `y2`.

    This function works exactly the same as matplotlib's fill_between, except
    that it also plots a proxy artist (specifically, a rectangle of 0 size)
    so that it can be added it appears on a legend.
    """
    ax = ax if ax is not None else plt.gca()
    ax.fill_between(x, y1, y2, interpolate=True, **kwargs)
    p = plt.Rectangle((0, 0), 0, 0, **kwargs)
    ax.add_patch(p)

# -----------------------------------------------------------------------------
def plotNumberFunction(sampleDirList=None, figDir="./", 
                       plotName="numberfunc", 
                       dataPortion="central"):

print "Plotting number function"

plt.clf()
plt.xlabel("$R_{eff}$ [$h^{-1}Mpc$]", fontsize=14)
plt.ylabel(r"log ($n$ (> R) [$h^3$ Gpc$^{-3}$])", fontsize=14)

for (iSample,sampleDir) in enumerate(sampleDirList):
  with open(workDir+sampleDir+"/sample_info.dat", 'rb') as input:
    sample = pickle.load(input)
  filename = sampleDir+"/centers_"+dataPortion+"_"+sample.fullName+".out"
  if not os.access(filename, os.F_OK):
    print "File not found: ", filename
  else:
    data = np.loadtxt(filename, comments="#")[:,4]

  if sample.dataType == "observation":
    # look for the mask file
    if os.access(sample.maskFile, os.F_OK):
      maskFile = sample.maskFile
    else:
      maskFile = sampleDir+"/"+os.path.basename(sample.maskFile)
      print "Using maskfile found in:", maskFile

    boxVol = vp.getSurveyProps(maskFile,
                               sample.zBoundary[0], sample.zBoundary[1],
                               sample.zRange[0], sample.zRange[1], "all",
                               selectionFuncFile=None)[0]
                               #selectionFuncFile=sample.selFunFile)[0]
    boxVol *= obsFudgeFactor
  else:
    boxVol = sample.boxLen*sample.boxLen*(sample.zBoundaryMpc[1] -
                                          sample.zBoundaryMpc[0])

  boxVol *= 1.e-9 # Mpc->Gpc

  bins = args.xmax/5.
  hist, binEdges = np.histogram(data, bins=bins, range=(0., 100.))
  binCenters = 0.5*(binEdges[1:] + binEdges[:-1])

  nvoids = len(data)
  var = hist * (1. - hist/nvoids)
  sig = np.sqrt(var)

  lowerbound = hist - sig
  upperbound = hist + sig

  mean = np.log10(hist/boxVol)
  lowerbound = np.log10(lowerbound/boxVol)
  upperbound = np.log10(upperbound/boxVol)

  lineColor = colorList[iSample]
  lineTitle = sample.fullName
  
  trim = (bounds[0] > .01)
  mean = mean[trim]
  binCentersToUse = binCenters[trim]
  lower = lowerbound[trim]
  upper = upperbound[trim]

  alpha = 0.55
  fill_between(binCentersToUse, lower, upper,
             label=lineTitle, color=lineColor,
             alpha=alpha,
             )

  lineStyle = '-'
  plt.plot(binCentersToUse, mean, lineStyle,
              color=lineColor,
              linewidth=3)

  plt.legend(loc = "upper right", fancybox=True, prop={'size':14})

  plt.savefig(figDir+"/fig_"+plotName+".pdf", bbox_inches="tight")
  plt.savefig(figDir+"/fig_"+plotName+".eps", bbox_inches="tight")
  plt.savefig(figDir+"/fig_"+plotName+".png", bbox_inches="tight")
