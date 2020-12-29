#!/usr/bin/env python
#+
#   VIDE -- Void IDentification and Examination -- ./python_tools/misc_util/figureOutMask.py
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

# build a low-resolution mask given the actual galaxy positions

import numpy as np
import healpy as hp
import os
import shutil
import glob
import sys

galFile = "path/to/galaxy/file"
outMaskFile = "path/to/output/mask"

nside = 128
npix = hp.nside2npix(nside)
mask = np.zeros((npix))

for line in open(galFile):
  line = line.split()
  RA  = np.float(line[3])
  Dec = np.float(line[4])
  z   = np.float(line[5])

  phi   = np.pi/180.*RA
  theta = Dec*np.pi/180.
  theta = np.pi/2. - Dec*np.pi/180.
  pos = np.zeros((3))

  pix = hp.ang2pix(nside, theta, phi)
  mask[pix] = 1.

print "Sky Fraction:", 1.0*len(mask[mask>0])/len(mask)

hp.write_map(outMaskFile, mask) 
