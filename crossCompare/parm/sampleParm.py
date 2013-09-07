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
