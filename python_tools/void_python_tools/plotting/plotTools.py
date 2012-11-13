__all__=['plotRedshiftDistribution', 'plotSizeDistribution', 'plot1dProfiles',
         'plotMarg1d', 'plotNumberDistribution']

from void_python_tools.backend.classes import *
from plotDefs import *
import numpy as np
import os
import pylab as plt
import void_python_tools.apTools as vp

# -----------------------------------------------------------------------------
def plotRedshiftDistribution(workDir=None, sampleList=None, figDir=None, 
                     plotNameBase="zdist", 
                     showPlot=False, dataPortion=None, setName=None):

  plt.clf()
  plt.xlabel("Redshift")
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
    nbins = np.ceil((zMax-zMin)/0.02)

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


# -----------------------------------------------------------------------------
def plotSizeDistribution(workDir=None, sampleList=None, figDir=None, 
                     plotNameBase="sizedist", 
                     showPlot=False, dataPortion=None, setName=None):

  plt.clf()
  plt.xlabel("Void Radius (Mpc/h)")
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

    xMin = 5
    xMax = 140

    range = (xMin, xMax)
    nbins = np.ceil((xMax-xMin)/5)

    #thisMax = np.max(data[:,5])
    #thisMin = np.min(data[:,5])
    #if thisMax > xMax: xMax = thisMax
    #if thisMin < xMin: xMin = thisMin

    plt.hist(data[:,4], bins=nbins,
             label=lineTitle, color=colorList[iSample],
             histtype = "step", range=range,
             linewidth=linewidth)

  plt.legend(title = "Samples", loc = "upper right")
  plt.title(plotTitle)

  plt.xlim(xMin, xMax)
  #plt.xlim(xMin, xMax*1.4) # make room for legend

  plt.savefig(figDir+"/fig_"+plotName+".pdf", bbox_inches="tight")
  plt.savefig(figDir+"/fig_"+plotName+".eps", bbox_inches="tight")
  plt.savefig(figDir+"/fig_"+plotName+".png", bbox_inches="tight")

  if showPlot:
    os.system("display %s" % figDir+"/fig_"+plotName+".png")


# -----------------------------------------------------------------------------
def plot1dProfiles(workDir=None, sampleList=None, figDir=None, 
                   plotNameBase="1dprofile", 
                   showPlot=False, dataPortion=None, setName=None):

  plt.clf()
  plt.xlabel(r"$R/R_{v,\mathrm{max}}$")
  plt.ylabel(r"$n / \bar n$")

  for (iSample,sample) in enumerate(sampleList):

    sampleName = sample.fullName

    for (iStack,stack) in enumerate(sample.stacks):

      plotTitle = setName
      plotName = plotNameBase

      runSuffix = getStackSuffix(stack.zMin, stack.zMax, stack.rMin, 
                                 stack.rMax, dataPortion)

      plotTitle = sampleName + ", z = "+str(stack.zMin)+"-"+str(stack.zMax)+", R = "+str(stack.rMin)+"-"+str(stack.rMax)+ r" $h^{-1}$ Mpc"

      filename = workDir+"/sample_"+sampleName+"/stacks_"+runSuffix+\
                 "profile_1d.txt"

      if not os.access(filename, os.F_OK):
        print "File not found: ", filename
        continue

      data = np.loadtxt(filename, comments="#")
      if data.ndim == 1:
        print " Too few!"
        continue

      data[:,1] /= stack.rMax
      plt.ylim(ymin=0.0, ymax=np.amax(data[:,2])+0.1)
      plt.xlim(xmin=0.0, xmax=2.1)
      plt.plot(data[:,1], data[:,2], label=lineTitle, color=colorList[0],
               linewidth=linewidth)

      plt.title(plotTitle)

      plt.savefig(figDir+"/fig_"+plotName+".pdf", bbox_inches="tight")
      plt.savefig(figDir+"/fig_"+plotName+".eps", bbox_inches="tight")
      plt.savefig(figDir+"/fig_"+plotName+".png", bbox_inches="tight")

      if showPlot:
        os.system("display %s" % figDir+"/fig_"+plotName+".png")


# -----------------------------------------------------------------------------
def plotMarg1d(workDir=None, sampleList=None, figDir=None, 
               plotNameBase="marg1d", 
               showPlot=False, dataPortion=None, setName=None):

    plotNames = ("Om", "w0", "wa")
    plotTitles = ("$\Omega_M$", "$w_0$", "$w_a$")
    files = ("Om", "w0", "wa")

    for iPlot in range(len(plotNames)):

      plt.clf()

      plotName = plotNameBase+"_"+plotNames[iPlot]+"_"+dataPortion
      plotTitle = plotTitles[iPlot]
      dataFile = workDir + "/likelihoods_"+dataPortion+"_"+files[iPlot]+".dat"

      plt.xlabel(plotTitle, fontsize="20")
      plt.ylabel("Likelihood", fontsize="20")
      plt.ylim(0.0, 1.0)

      data = np.loadtxt(dataFile, comments="#")
      plt.plot(data[:,0], data[:,1], color='k', linewidth=linewidth)

      plt.savefig(figDir+"/fig_"+plotName+".pdf", bbox_inches="tight")
      plt.savefig(figDir+"/fig_"+plotName+".eps", bbox_inches="tight")
      plt.savefig(figDir+"/fig_"+plotName+".png", bbox_inches="tight")

      if showPlot:
        os.system("display %s" % figDir+"/fig_"+plotName+".png")



# -----------------------------------------------------------------------------
def plotNumberDistribution(workDir=None, sampleList=None, figDir=None, 
                     plotNameBase="numberdist", 
                     showPlot=False, dataPortion=None, setName=None):

  plt.clf()
  plt.xlabel("Void Radius (Mpc/h)")
  plt.ylabel(r"N > R [h^3 Mpc^{-3}]")

  plotTitle = setName

  plotName = plotNameBase

  plt.yscale('log')

  for (iSample,sample) in enumerate(sampleList):

    sampleName = sample.fullName
    lineTitle = sampleName

    if sample.dataType == "observation":
      boxVol = vp.getSurveyProps(sample.maskFile, stack.zMin, stack.zMax, 
                                 sample.zRange[0], sample.zRange[1], "all",
                                 selectionFuncFile=sample.selFunFile)[0]
    else:
      boxVol = sample.boxLen*sample.boxLen*(sample.zBoundaryMpc[1] - 
                                            sample.zBoundaryMpc[0])

    filename = workDir+"/sample_"+sampleName+"/centers_"+dataPortion+"_"+\
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

  plt.legend(title = "Samples", loc = "upper right")
  plt.title(plotTitle)

  plt.xlim(xMin, xMax)
  #plt.xlim(xMin, xMax*1.4) # make room for legend

  plt.savefig(figDir+"/fig_"+plotName+".pdf", bbox_inches="tight")
  plt.savefig(figDir+"/fig_"+plotName+".eps", bbox_inches="tight")
  plt.savefig(figDir+"/fig_"+plotName+".png", bbox_inches="tight")

  if showPlot:
    os.system("display %s" % figDir+"/fig_"+plotName+".png")



