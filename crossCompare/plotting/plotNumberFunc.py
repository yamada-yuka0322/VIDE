#+
#   VIDE -- Void IDEntification pipeline -- ./crossCompare/plotting/plotNumberFunc.py
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

# plots cumulative distributions of number counts

from void_python_tools.backend import *
from void_python_tools.plotting import *
import void_python_tools.apTools as vp
import imp
import pickle
import os
import matplotlib.pyplot as plt
import numpy as np
import argparse
from scipy.stats import ks_2samp

# ------------------------------------------------------------------------------

plotNameBase = "compdist"

obsFudgeFactor = 1.0 # what fraction of the volume are we *reall* capturing?
#obsFudgeFactor = .66 # what fraction of the volume are we *reall* capturing?

linewidth = 1

parser = argparse.ArgumentParser(description='Plot.')
parser.add_argument('--show', dest='showPlot', action='store_const',
                   const=True, default=False,
                   help='display the plot (default: just write eps)')
parser.add_argument('--parmFile', dest='parmFile', default='datasetsToPlot.py',
                    help='path to parameter file')
args = parser.parse_args()

nErrorBars = 10
plotMax = 120
errorBarsX = np.linspace(0, plotMax, num=nErrorBars)

# ------------------------------------------------------------------------------

filename = args.parmFile
print " Loading parameters from", filename
if not os.access(filename, os.F_OK):
  print "  Cannot find parameter file %s!" % filename
  exit(-1)
parms = imp.load_source("name", filename)
globals().update(vars(parms))

if not os.access(figDir, os.F_OK):
  os.makedirs(figDir)

dataSampleList = []

for sampleDir in sampleDirList:
  with open(workDir+sampleDir+"/sample_info.dat", 'rb') as input:
    dataSampleList.append(pickle.load(input))

plt.clf()
plt.xlabel("Void Radius [Mpc/h]")
plt.ylabel(r"N > R [$h^3$ Gpc$^{-3}$]")
plt.yscale('log')
plt.xlim(xmax=plotMax)

plotName = plotNameBase
allData = []

for dataPortion in dataPortions: 
  print "Data portion:", dataPortion
  sizeList = []
  for (iSample,sample) in enumerate(dataSampleList):

    sampleName = sample.fullName
    lineTitle = sampleName

    if sample.dataType == "observation":
      boxVol = vp.getSurveyProps(sample.maskFile, 
                                 sample.zBoundary[0], sample.zBoundary[1], 
                                 sample.zRange[0], sample.zRange[1], "all",
                                 selectionFuncFile=None)[0]
                                 #selectionFuncFile=sample.selFunFile)[0]
      boxVol *= obsFudgeFactor
    else:
      boxVol = sample.boxLen*sample.boxLen*(sample.zBoundaryMpc[1] - 
                                          sample.zBoundaryMpc[0])

    boxVol *= 1.e-9 # Mpc->Gpc

    filename = workDir+"/"+sampleDirList[iSample]+"/centers_"+dataPortion+"_"+\
               sampleName+".out"
    if not os.access(filename, os.F_OK):
      print "File not found: ", filename
      continue

    data = np.loadtxt(filename, comments="#")
    if data.ndim == 1:
      print " Too few!"
      continue
    data = data[:,4]
    indices = np.arange(0, len(data), 1)
    numVoids = indices[::-1]/boxVol
    voidSizes = np.sort(data)

    sizeList.append(voidSizes)
    myErrorBarsX = []
    errorBarsY = []
    errorBarsDY = []
    errorBarsIdx = []
    for errorBarLoc in errorBarsX:
      nearestIdx = (np.abs(voidSizes-errorBarLoc)).argmin()
      if nearestIdx == 0: continue
      myErrorBarsX.append(errorBarLoc)
      errorBarsIdx.append(nearestIdx)
      errorBarsY.append(numVoids[nearestIdx])
      errorBarsDY.append(np.sqrt(numVoids[nearestIdx]))

    thisPlot = plt.plot(voidSizes, numVoids, '-',
             color=colorList[iSample],
             linewidth=linewidth, label=lineTitle)
    plt.errorbar(myErrorBarsX, errorBarsY, errorBarsDY,
             ecolor=colorList[iSample],
             fmt=None, label='_nolegend_', capsize=0)

    hist, bin_edges = np.histogram(data, bins=100, range=(0,100)) 
    allData.append(hist)

  plt.legend(title = "Samples", loc = "upper right", prop={'size':8})
  #plt.title(plotTitle)

  # compute K-S statistic for all pairs of sets
  for (i,sample1) in enumerate(dataSampleList):
    for (j,sample2) in enumerate(dataSampleList):
      if j <= i: continue
      ks, pval = ks_2samp(sizeList[i][:], sizeList[j][:])
      print sample1.fullName, sample2.fullName, pval


plt.savefig(figDir+"/fig_"+plotName+".pdf", bbox_inches="tight")
plt.savefig(figDir+"/fig_"+plotName+".eps", bbox_inches="tight")
plt.savefig(figDir+"/fig_"+plotName+".png", bbox_inches="tight")

dataFile = figDir+"/data_"+plotName+".dat"
fp = open(dataFile, 'w')
fp.write("# R [Mpc/h], N [h^3 Gpc^-3]\n")
fp.write("# ")
for sample in dataSampleList:
  fp.write(sample.fullName+" ")
fp.write("\n")
for i in xrange(100):
  fp.write(str(bin_edges[i]) + " ")
  for iSample in xrange(len(dataSampleList)):
    fp.write(str(allData[iSample][i])+" ")
  fp.write("\n")
fp.close()

if args.showPlot:
  os.system("display %s" % figDir+"/fig_"+plotName+".png")


