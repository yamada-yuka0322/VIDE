#+
#   VIDE -- Void IDEntification pipeline -- ./crossCompare/analysis/datasetsToAnalyze.py
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


workDir = "/home/psutter2/workspace/Voids/"
dataDir = "/home/psutter2/workspace/Voids/crossCompare/mergerTree/"
 

CTOOLS_PATH = "../../c_tools/"

baseSampleDir = "mergertree512/mt_ss0.1/sample_md_ss0.1_z0.00_d00/"

sampleDirList = [ 
                  "mergertree512/mt_ss1e-05/sample_md_ss1e-05_z0.00_d00/",
                  #"mergertree512/mt_ss0.000175/sample_md_ss0.000175_z0.00_d00/",
                  #"mergertree512/mt_ss0.0004/sample_md_ss0.0004_z0.00_d00/",
                  #"mergertree512/mt_ss0.001/sample_md_ss0.001_z0.00_d00/",
                  #"mergertree512/mt_ss0.002/sample_md_ss0.002_z0.00_d00/",
                  "mergertree512/mt_ss0.01/sample_md_ss0.01_z0.00_d00/",
                  "mergertree512/mt_hod_dr72dim2/sample_md_hod_dr72dim2_z0.00_d00/",
                  "mergertree512/mt_hod_dr9mid/sample_md_hod_dr9mid_z0.00_d00/",
                  "mergertree512/mt_halos_min1.2e+13/sample_md_halos_min1.2e+13_z0.00_d00/",
                ]

dataPortion = "central"

