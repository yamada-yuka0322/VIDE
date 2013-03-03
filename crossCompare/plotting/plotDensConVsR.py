#+
#   VIDE -- Void IDEntification pipeline -- ./crossCompare/plotting/plotDensConVsR.py
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

# plots cumulative distributions of number counts versus density contrast

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

plotNameBase = "densconvsr"

obsFudgeFactor = .66 # what fraction of the volume are we *reall* capturing?

parser = argparse.ArgumentParser(description='Plot.')
parser.add_argument('--show', dest='showPlot', action='store_const',
                   const=True, default=False,
                   help='display the plot (default: just write eps)')
parser.add_argument('--parmFile', dest='parmFile', default='datasetsToPlot.py',
                    help='path to parameter file')
args = parser.parse_args()

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
plt.ylabel("Void Density Contrast")
plt.xlabel("Void Radius [Mpc/h]")
#plt.xlim(xmax=5.)
plt.yscale('log')

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

  print " Loading", sampleName
  filename = workDir+"/"+sampleDirList[iSample]+"/centers_"+dataPortion+"_"+\
             sampleName+".out"
  if not os.access(filename, os.F_OK):
    print "File not found: ", filename
    continue

  data = np.loadtxt(filename, comments="#")
  if data.ndim == 1:
    print " Too few!"
    continue

  plt.plot(data[:,4], data[:,8], '-',
           label=lineTitle, color=colorList[iSample],
           linewidth=linewidth)

plt.legend(title = "Samples", loc = "upper right", prop={'size':8})
#plt.title(plotTitle)

plt.savefig(figDir+"/fig_"+plotName+".pdf", bbox_inches="tight")
plt.savefig(figDir+"/fig_"+plotName+".eps", bbox_inches="tight")
plt.savefig(figDir+"/fig_"+plotName+".png", bbox_inches="tight")

if args.showPlot:
  os.system("display %s" % figDir+"/fig_"+plotName+".png")


