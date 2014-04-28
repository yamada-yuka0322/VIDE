#+
#   VIDE -- Void IDentification and Examination -- ./analysis/datasetsToAnalyze.py
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


outputDir = "/home/psutter2/workspace/Voids/analysis/xcor/"

# Sim parameters
Ni = 1000237  # Number of dark matter particles per dimension in simulation
ss = 0.1  # Subsampling fraction of dark matter particles to read
Mpart = 8.721e9  # Particle mass [M_sol]
Npart = int(ss*Ni)  # Particle number of subsample
Lboxcut = 0.  # Size of optional margin to be cut from the box [h^(-1)Mpc]
Nmesh = 256  # Interpolation meshlength
Nbin = 70  # Number of bins for power spectrum and correlation function
r_H = 3000.  # Hubble scale [h^(-1)Mpc]
ns = 0.95  # Spectral index
sigma_8 = 0.82  # Sigma_8
h = 0.7  # Dimensionless Hubble parameter

# Input files
matterDir = '/home/psutter2/workspace/Voids/catalogs/mergertree512/'
haloDir = '/home/psutter2/workspace/Voids/catalogs/mergertree512/'
matterFilename = 'mf_4s_1G_512_1.000'
haloFilename = 'mf_4s_1G_512_bgc2_1.000.sdf'
voidBaseDir = "/home/psutter2/workspace/Voids/"
 

sampleDirList = [ 
                  "mergertree512/mt_ss0.01/sample_mt_ss0.01_z0.00_d00/",
                ]

dataPortion = "central"

