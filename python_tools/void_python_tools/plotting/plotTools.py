__all__=['plotNumberCounts']

from void_python_tools.backend.classes import *
from plotDefs import *
import numpy as np
import os
import pylab as plt

# -----------------------------------------------------------------------------
def plotNumberCounts(workDir=None, sampleList=None, figDir=None, plotNameBase="numbercount", 
                     showPlot=False, dataPortion=None, setName=None):

  plt.clf()
  plt.xlabel("Comoving Distance (Mpc/h)")
  plt.ylabel("Number of Voids")

  plotTitle = setName

  plotName = plotNameBase

  xMin = 1.e00
  xMax = 0

  for (iSample,sample) in enumerate(sampleList):

    sampleName = sample.fullName
    lineTitle = sampleName

    filename = workDir+"/sample_"+sampleName+"/centers_"+dataPortion+"_"+\
               sampleName+".out"
    if not os.access(filename, os.F_OK):
      print "File not found: ", filename
      continue

    data = np.loadtxt(filename, comments="#")
    if data.ndim == 1:
      print " Too few!"
      continue

    zMin = sample.zRange[0]
    zMax = sample.zRange[1]

    range = (zMin, zMax)
    nbins = np.ceil((zMax-zMin)/0.1)

    thisMax = np.max(data[:,5])
    thisMin = np.min(data[:,5])
    if thisMax > xMax: xMax = thisMax
    if thisMin < xMin: xMin = thisMin

    plt.hist(data[:,5], bins=nbins,
             label=lineTitle, color=colorList[iSample],
             histtype = "step", range=range,
             linewidth=linewidth)

  #plt.legend(title = "Samples", loc = "upper right")
  plt.title(plotTitle)

  plt.xlim(xMin, xMax)
  #plt.xlim(xMin, xMax*1.4) # make room for legend

  plt.savefig(figDir+"/fig_"+plotName+".pdf", bbox_inches="tight")
  plt.savefig(figDir+"/fig_"+plotName+".eps", bbox_inches="tight")
  plt.savefig(figDir+"/fig_"+plotName+".png", bbox_inches="tight")

  if showPlot:
    os.system("display %s" % figDir+"/fig_"+plotName+".png")


