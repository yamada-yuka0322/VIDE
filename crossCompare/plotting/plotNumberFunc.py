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

# ------------------------------------------------------------------------------

obsFudgeFactor = 1.0 # what fraction of the volume are we *reall* capturing?
#obsFudgeFactor = .66 # what fraction of the volume are we *reall* capturing?

histBinWidth = 1 # Mpc

parser = argparse.ArgumentParser(description='Plot.')
parser.add_argument('--show', dest='showPlot', action='store_const',
                   const=True, default=False,
                   help='display the plot (default: just write eps)')
parser.add_argument('--binned', dest='binned', action='store_const',
                   const=True, default=False,
                   help='plot binned function (default: cumulative)')
parser.add_argument('--parm', dest='parm', default='datasetsToPlot.py',
                    help='path to parameter file')
parser.add_argument('--xmax', dest='xmax', default=120.,
                    help='x limit of plot')
parser.add_argument('--xmin', dest='xmin', default=20.,
                    help='x limit of plot')
args = parser.parse_args()

# ------------------------------------------------------------------------------


def loadData(sampleDir, dataPortion):
  with open(workDir+sampleDir+"/sample_info.dat", 'rb') as input:
    sample = pickle.load(input)

  filename = workDir+"/"+sampleDir+"/centers_"+dataPortion+"_"+sample.fullName+".out"
  if not os.access(filename, os.F_OK):
    print "File not found: ", filename
    return -1, -1, -1

  data = np.loadtxt(filename, comments="#")
  if data.ndim == 1:
    print " Too few!"
    return -1, -1, -1

  data = data[:,4]
  indices = np.arange(0, len(data), 1)
  sorted = np.sort(data)

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

  indices /= boxVol

  #xmin = sorted[0]
  #xmax = sorted[-1]
  #bins = int((xmax-xmin)/histBinWidth)
  bins = args.xmax/histBinWidth
  hist, binEdges = np.histogram(sorted, bins=bins, range=(0., args.xmax)) 
  #hist, binEdges = np.histogram(sorted, bins=bins, range=(xmin,xmax)) 
  binCenters = 0.5*(binEdges[1:] + binEdges[:-1])

  if not args.binned:
    foundStart = False
    for iBin in xrange(len(hist)):
      if not foundStart and hist[iBin] == 0:
        continue
      foundStart = True
      hist[iBin] = np.sum(hist[iBin:])

  hist /= boxVol

  hist = np.log10(hist)

  lineTitle = sample.nickName

  return hist, binCenters, lineTitle

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

# ------------------------------------------------------------------------------

print "Plotting number function"

filename = args.parm
print " Loading parameters from", filename
if not os.access(filename, os.F_OK):
  print "  Cannot find parameter file %s!" % filename
  exit(-1)
parms = imp.load_source("name", filename)
globals().update(vars(parms))

if not os.access(figDir, os.F_OK):
  os.makedirs(figDir)

plt.clf()
plt.xlabel("Void Radius [Mpc/h]")
plt.ylabel(r"log (N > R [$h^3$ Gpc$^{-3}$])")
#plt.yscale('log')
#plt.xlim(xmax=args.xmax)
#plt.xlim(xmin=args.xmin)
#plt.ylim(ymin=-1)
#plt.ylim(ymax=6)

plotNameBase = "numberfunc"
plotName = plotNameBase + "_" + plotLabel

for (iSample,sampleDir) in enumerate(sampleDirList):
 for dataPortion in dataPortions:
  # get all the data
  allHist = []
  if "ZZZZ" in sampleDir:
    for fileZ in fileList:
      thisSampleDir = sampleDir.replace("ZZZZ", fileZ)
      hist, binCenters, lineTitle = loadData(thisSampleDir, dataPortion)
      if lineTitle == -1: continue
      allHist.append(hist)
   
    lineLabel = lineTitle.replace(fileZ, "all")
    lineLabel += ", " + dataPortion
    
    maxHist = 1.*allHist[-1]
    minHist = 1.*allHist[-1]
    for iHist in xrange(len(allHist)-1):
      maxHist = np.maximum(maxHist, allHist[iHist])
      minHist = np.minimum(minHist, allHist[iHist])

    trim = (maxHist > 1) 
    minHist = minHist[trim]
    maxHist = maxHist[trim]
    binCentersToUse = binCenters[trim]
    if dataPortion == "central":
      hatch = '/'
    else:
      hatch = None
    fill_between(binCentersToUse, minHist, maxHist,
             label=lineLabel, color=colorList[iSample],
             alpha=0.5, hatch=hatch
             )

  else:

    hist, binCenters, lineLabel = loadData(sampleDir, dataPortion)
    trim = (hist > 1) 
    hist = hist[trim]
    binCentersToUse = binCenters[trim]
    if lineLabel == -1: continue
    lineLabel += ", " + dataPortion
    if dataPortion == "central":
      lineStyle = '--'
    else:
      lineStyle = '-'
    plt.plot(binCentersToUse, hist, lineStyle,
              label=lineLabel, color=colorList[iSample],
              linewidth=linewidth)
    

plt.legend(title = "Samples", loc = "upper right", prop={'size':8})
plt.title("Number func - "+plotTitle)

plt.savefig(figDir+"/fig_"+plotName+".pdf", bbox_inches="tight")
plt.savefig(figDir+"/fig_"+plotName+".eps", bbox_inches="tight")
plt.savefig(figDir+"/fig_"+plotName+".png", bbox_inches="tight")

if args.showPlot:
  os.system("display %s" % figDir+"/fig_"+plotName+".png")


