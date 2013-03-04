#!/usr/bin/env python


workDir = "/home/psutter2/workspace/Voids/"
dataDir = "/home/psutter2/workspace/Voids/crossCompare/mergerTree/"
 

CTOOLS_PATH = "../../c_tools/"

baseSampleDir = "multidark/md_ss0.05/sample_md_ss0.05_z0.56_d00/"

sampleDirList = [ 
                  "multidark/md_ss1e-05/sample_md_ss1e-05_z0.56_d00/",
                  #"multidark/md_ss0.000175/sample_md_ss0.000175_z0.56_d00/",
                  #"multidark/md_ss0.0004/sample_md_ss0.0004_z0.56_d00/",
                  #"multidark/md_ss0.001/sample_md_ss0.001_z0.56_d00/",
                  #"multidark/md_ss0.002/sample_md_ss0.002_z0.56_d00/",
                  "multidark/md_ss0.01/sample_md_ss0.01_z0.56_d00/",
                  "multidark/md_hod_dr72dim2/sample_md_hod_dr72dim2_z0.56_d00/",
                  "multidark/md_hod_dr9mid/sample_md_hod_dr9mid_z0.56_d00/",
                  "multidark/md_halos_min1.2e+13/sample_md_halos_min1.2e+13_z0.56_d00/",
                ]

dataPortion = "central"

