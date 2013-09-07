#!/usr/bin/env python

# plots cumulative distributions of number counts

import matplotlib
matplotlib.use('Agg')
from void_python_tools.backend import *
from void_python_tools.plotting import *
import void_python_tools.apTools as vp
import imp
import pickle
import os
import matplotlib.pyplot as plt
import numpy as np
import argparse
from globalOptions import *
from scipy.optimize import curve_fit

# plots the ratio of sizes for matched voids

# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='Plot.')
parser.add_argument('--show', dest='showPlot', action='store_const',
                   const=True, default=False,
                   help='display the plot (default: just write eps)')
parser.add_argument('--parm', dest='parm', default='datasetsToPlot.py',
                    help='path to parameter file')
args = parser.parse_args()

# -----------------------------------------------------------------------------
def Linear(x, a, b):
  return a*x+b

# ------------------------------------------------------------------------------

print "Plotting match size ratio"

filename = args.parm
print " Loading parameters from", filename
if not os.access(filename, os.F_OK):
  print "  Cannot find parameter file %s!" % filename
  exit(-1)
parms = imp.load_source("name", filename)
globals().update(vars(parms))

if not os.access(figDir, os.F_OK):
  os.makedirs(figDir)

dataSampleList = []
compareSampleList = []

with open(workDir+baseSampleDir+"/sample_info.dat", 'rb') as input:
  baseSample = pickle.load(input)

for sampleDir in sampleDirList:
  with open(workDir+sampleDir+"/sample_info.dat", 'rb') as input:
    thisSample = pickle.load(input)
    if compareSampleTag in thisSample.fullName:
      compareSampleList.append(thisSample)
    else:
      dataSampleList.append(thisSample)

plt.clf()
#plt.yscale('log')

numSubPlots = len(dataSampleList)
fig, axesList = plt.subplots(numSubPlots, sharex=True, sharey=True)
axesList = np.atleast_1d(axesList)

for (iSample,sample) in enumerate(dataSampleList):

  if sample.fullName == baseSample.fullName: continue

  sampleName = sample.fullName
  lineTitle = sample.nickName[:-10]

  plotNameBase = "matchvolrelradius"
  plotName = plotNameBase + "_" + plotLabel# + "_" + sampleName

  filename = outputDir+"/voidOverlap_"+baseSample.fullName+"_"+sampleName+"_summary.out"
  if not os.access(filename, os.F_OK):
    print "File not found: ", filename
    continue

  data = np.loadtxt(filename, comments="#")
  if data.ndim == 1:
    print " Too few!"
    continue

  # find a sample to compare it to
  for compareSample in compareSampleList:
    if compareSample.nickName[:10] == sample.nickName[:10]:
      filename = outputDir+"/voidOverlap_"+baseSample.fullName+"_"+compareSample.fullName+"_summary.out"
      compareData = np.loadtxt(filename, comments="#")
      axesList[iSample].scatter(compareData[:,1], compareData[:,2],
            color='blue', alpha=alpha, s=pointsize)
 
  axesList[iSample].scatter(data[:,1], data[:,2],
           label=lineTitle, color='red', alpha=alpha, s=pointsize)

  try:
    popt, pcov = curve_fit(Linear, data[:,1],
                           data[:,2], maxfev=10000, xtol=5.e-3)
  except RuntimeError:
    print "Warning: no convergence reached"

  label = r"$a$=%.2f, $b$=%.2f" % (popt[0], popt[1])
  radii = np.arange(rMin, rMax, 1)
  relRad = Linear(radii, popt[0], popt[1])
  #axesList[iSample].plot(radii, relRad, '-', label=label,
  #         color='black', linewidth=2)

  axesList[iSample].legend(loc = "best", fancybox=True,prop={'size':10})

  plt.ylim(0., 2.5)
     
  yticks = axesList[iSample].yaxis.get_major_ticks()
  yticks[-1].label1.set_visible(False)
  yticks[0].label1.set_visible(False)
 
  axesList[iSample].set_xlim([rMin,rMax])


fig.subplots_adjust(hspace=0)

plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

axesList[0].set_title(plotTitle, fontsize=14)
#axesList[0].set_title("Match size ratio - "+plotTitle)
fig.text(0.5, 0.04, r'$R_{eff}$ [$h^{-1}$Mpc]', ha='center', va='center', fontsize=14)
fig.text(0.06, 0.5, 'Match Relative Radius', ha='center', va='center', rotation='vertical', fontsize=14)

#plt.legend(title = "Samples", loc = "upper left", prop={'size':8})
#plt.title("Match size ratio - "+plotTitle+" - "+sample.nickName)

plt.savefig(figDir+"/fig_"+plotName+".pdf", bbox_inches="tight")
plt.savefig(figDir+"/fig_"+plotName+".eps", bbox_inches="tight")
plt.savefig(figDir+"/fig_"+plotName+".png", bbox_inches="tight")

if args.showPlot:
  os.system("display %s" % figDir+"/fig_"+plotName+".png")


