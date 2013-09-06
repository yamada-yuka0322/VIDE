#+
#   VIDE -- Void IDEntification pipeline -- ./python_tools/void_python_tools/plotting/plotTools.py
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
import void_python_tools.xcor as xcor

# -----------------------------------------------------------------------------
def plotRedshiftDistribution(workDir=None, sampleList=None, figDir=None, 
                     plotNameBase="zdist", 
                     showPlot=False, dataPortion=None, setName=None):

  plt.ioff()
  maptplotlib.pyplot.switch_backed('Agg')

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


  plt.ioff()
  maptplotlib.pyplot.switch_backed('Agg')

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


  plt.ioff()
  maptplotlib.pyplot.switch_backed('Agg')

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

    plt.ioff()
    maptplotlib.pyplot.switch_backed('Agg')

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

  plt.ioff()
  maptplotlib.pyplot.switch_backed('Agg')

  plt.clf()
  plt.xlabel("Void Radius (Mpc/h)")
  plt.ylabel(r"N > R [$h^3$ Mpc$^{-3}$]")

  plotTitle = setName

  plotName = plotNameBase

  plt.yscale('log')

  for (iSample,sample) in enumerate(sampleList):

    sampleName = sample.fullName
    lineTitle = sampleName

    if sample.dataType == "observation":
      boxVol = vp.getSurveyProps(sample.maskFile, 
                                 sample.zBoundary[0], sample.zBoundary[1],
                                 sample.zRange[0], sample.zRange[1], "all",
                                 selectionFuncFile=None)[0]
    else:
      boxVol = sample.boxLen*sample.boxLen*(sample.zBoundaryMpc[1] - 
                                            sample.zBoundaryMpc[0])

    boxVol *= 1.e-9

    filename = workDir+"/sample_"+sampleName+"/untrimmed_centers_"+dataPortion+"_"+\
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

  plt.savefig(figDir+"/fig_"+plotName+".pdf", bbox_inches="tight")
  plt.savefig(figDir+"/fig_"+plotName+".eps", bbox_inches="tight")
  plt.savefig(figDir+"/fig_"+plotName+".png", bbox_inches="tight")

  if showPlot:
    os.system("display %s" % figDir+"/fig_"+plotName+".png")

# -----------------------------------------------------------------------------
def plotVoidDistribution(workDir=None, sampleList=None, figDir=None, 
                     plotNameBase="dv", 
                     showPlot=False, dataPortion=None, setName=None):

  plt.ioff()
  maptplotlib.pyplot.switch_backed('Agg')
  Nmesh = 256
 
  for (iSample,sample) in enumerate(sampleList):
    plt.clf()
    plt.xlabel("x (Mpc/h)")
    plt.ylabel("y (Mpc/h)")

    sampleName = sample.fullName
  
    plotName = plotNameBase+"_"+sampleName

    filename = workDir+"/sample_"+sampleName+"/untrimmed_centers_"+dataPortion+"_"+\
               sampleName+".out"

    if not os.access(filename, os.F_OK):
      print "File not found: ", filename
      continue

    void_file = open(filename,'r')
    void_header1 = void_file.readline().split(",")
    void_data1 = np.reshape(void_file.read().split(),
                            (-1,len(void_header1))).astype(np.float32)
    void_file.close()

    xv = void_data1[:,0:3]
    rv = void_data1[:,4]
    vv = void_data1[:,6]
    dv, wm, ws = xcor.cic(xv, sample.boxLen, Lboxcut = 0., Nmesh = Nmesh)
    plt.imshow(np.sum(dv+1,2)/Nmesh,extent=[0,sample.boxLen,0,sample.boxLen],
               aspect='equal',
               cmap='YlGnBu_r',interpolation='gaussian')
    plt.colorbar()

    plt.title("Voids: "+sampleName)

    plt.savefig(figDir+"/fig_"+plotName+".pdf", bbox_inches="tight")
    plt.savefig(figDir+"/fig_"+plotName+".eps", bbox_inches="tight")
    plt.savefig(figDir+"/fig_"+plotName+".png", bbox_inches="tight")

    if showPlot:
      os.system("display %s" % figDir+"/fig_"+plotName+".png")

