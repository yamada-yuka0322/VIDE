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
__all__=['plotNumberFunction',]

from void_python_tools.backend.classes import *
from plotDefs import *
import numpy as np
import os
import pylab as plt
import void_python_tools.apTools as vp
from void_python_tools.voidUtil import getArray

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
def plotNumberFunction(catalogList, 
                       figDir="./", 
                       plotName="numberfunc",
                       cumulative=True,
                       binWidth=1):

# plots a cumulative number function
#   catalogList: list of void catalogs to plot
#   figDir: output directory for figures
#   plotName: name to prefix to all outputs
#   cumulative: if True, plots cumulative number function
#   binWidth: width of histogram bins in Mpc/h
# returns:
#   numberFuncList: array of len(catalogList), 
#                   each element has array of size bins, number, +/- 1 sigma

  print "Plotting number function"
  
  plt.clf()
  plt.xlabel("$R_{eff}$ [$h^{-1}Mpc$]", fontsize=14)
  
  if cumulative:
    plt.ylabel(r"log ($n$ (> R) [$h^3$ Gpc$^{-3}$])", fontsize=14)
  else:
    plt.ylabel(r"log ($dn/dR$ [$h^3$ Gpc$^{-3}$])", fontsize=14)
 
  numberFuncList = [] 

  for (iSample,catalog) in enumerate(catalogList):
    sample = catalog.sampleInfo
    data = getArray(catalog.voids, 'radius')
  
    if sample.dataType == "observation":
      maskFile = sample.maskFile
  
      boxVol = vp.getSurveyProps(maskFile,
                                 sample.zBoundary[0], sample.zBoundary[1],
                                 sample.zRange[0], sample.zRange[1], "all",
                                 selectionFuncFile=sample.selFunFile)[0]
    else:
      boxVol = sample.boxLen*sample.boxLen*(sample.zBoundaryMpc[1] -
                                            sample.zBoundaryMpc[0])
  
    boxVol *= 1.e-9 # Mpc->Gpc
  
    bins = 100./binWidth
    hist, binEdges = np.histogram(data, bins=bins, range=(0., 100.))
    binCenters = 0.5*(binEdges[1:] + binEdges[:-1])

    if cumulative:
      foundStart = False
      for iBin in xrange(len(hist)):
        if not foundStart and hist[iBin] == 0:
          continue
        foundStart = True
        hist[iBin] = np.sum(hist[iBin:])
  
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
    
    trim = (lowerbound > .01)
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

    numberFuncList.append((binCentersToUse, mean, lower, upper))

  return numberFuncList 
