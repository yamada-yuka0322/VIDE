#+
#   VIDE -- Void IDentification and Examination -- ./crossCompare/parm/sampleParm.py
#   Copyright (C) 2010-2014 Guilhem Lavaux
#   Copyright (C) 2011-2014 P. M. Sutter
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


workDir = "" # base directory for all samples
outputDir = "" 
logDir  = "./logs/"
figDir  = "./figs/"

# path to c_tools directory in VIDE
CTOOLS_PATH = "/home/psutter2/projects/Voids/vide/c_tools/"

# the path under workDir/ which which holds the sample you want to compare againt (e.g., the fiducial case)
baseSampleDir = "sample_/"

# comma-separated list of samples to compare against baseSampleDir
sampleDirList = [ 
                  "sample1_/",
                  "sample2_/",
                  "sample3_/",
                ]

dataPortions = [ "all" ]

# this name gets appended to all output filenames and plots
plotLabel = "test"

# title to place on plots
plotTitle = "Test"

# don't touch this for now; it will be fully implemented and explained later
compareSampleTag = ""
doTheory = False
