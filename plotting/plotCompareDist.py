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

from datasetsToPlot import *

plotNameBase = "compdist"

obsFudgeFactor = .66 # what fraction of the volume are we *reall* capturing?

parser = argparse.ArgumentParser(description='Plot.')
parser.add_argument('--show', dest='showPlot', action='store_const',
                   const=True, default=False,
                   help='display the plot (default: just write eps)')
args = parser.parse_args()

# ------------------------------------------------------------------------------

if not os.access(figDir, os.F_OK):
  os.makedirs(figDir)

dataSampleList = []

for sampleDir in sampleDirList:
  with open(workDir+sampleDir+"/sample_info.dat", 'rb') as input:
    dataSampleList.append(pickle.load(input))

plt.clf()
plt.xlabel("Void Radius (Mpc/h)")
plt.ylabel(r"N > R [$h^3$ Gpc$^{-3}$]")
plt.yscale('log')
plt.xlim(xmax=80.)

plotName = plotNameBase

for (iSample,sample) in enumerate(dataSampleList):

  sampleName = sample.fullName
  lineTitle = sampleName

  if sample.dataType == "observation":
    boxVol = vp.getSurveyProps(sample.maskFile, 
                               sample.zBoundary[0], sample.zBoundary[1], 
                               sample.zRange[0], sample.zRange[1], "all",
                               selectionFuncFile=sample.selFunFile)[0]
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
  sorted = np.sort(data)

  plt.plot(sorted, indices[::-1]/boxVol, '-',
           label=lineTitle, color=colorList[iSample],
           linewidth=linewidth)

plt.legend(title = "Samples", loc = "upper right", prop={'size':8})
#plt.title(plotTitle)

plt.savefig(figDir+"/fig_"+plotName+".pdf", bbox_inches="tight")
plt.savefig(figDir+"/fig_"+plotName+".eps", bbox_inches="tight")
plt.savefig(figDir+"/fig_"+plotName+".png", bbox_inches="tight")

if args.showPlot:
  os.system("display %s" % figDir+"/fig_"+plotName+".png")


