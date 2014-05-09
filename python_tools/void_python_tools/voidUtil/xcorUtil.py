#+
#   VIDE -- Void IDentification and Examination -- ./analysis/xcor.py
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
#+
#   VIDE -- Void IDentification and Examination -- ./crossCompare/analysis/mergerTree.py
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

from void_python_tools.backend import *
import imp
import pickle
import argparse
import os
import string
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
from matplotlib.ticker import NullFormatter
import random
import sys
__all__=['computeCrossCor',]

# ------------------------------------------------------------------------------

def computeCrossCor(catalogDir, 
                    outputDir="./", logDir="./",
                    matchMethod="useID", dataPortion="central",
                    strictMatch=True,
                    pathToCTools="../../../c_tools"):

# Computes void-void and void-matter(galaxy) correlations
#  baseCatalogDir: directory of catalog 
#  compareCatagDir: directory of catalog 2
#  outputDir: directory to place outputs
#  logDir: directory to place log files
#  matchMethod: "useID" to use unique IDs, "prox" to use overlap of Voronoi cells
#  dataPortion: "central" or "all"
#  strictMatch: if True, only attempt to match to trimmed catalog
#  pathToCTools: path to location of VIDE c_tools directory

  if not os.access(outputDir, os.F_OK):
    os.makedirs(outputDir)

  with open(catalogDir+"/sample_info.dat", 'rb') as input:
    sample = pickle.load(input)

  print " Working with", sample.fullName, "...",
  sys.stdout.flush()

  sampleName = sample.fullName

  # Sim parameters
  Lbox = sample.boxLen  # Boxlength [h^(-1)Mpc]
  Lbox -= 2*Lboxcut  # Reduced boxlength [h^(-1)Mpc]
  Om = sample.omegaM  # Omega_m
  Ol = 1.-Om  # Omega_l
  z = sample.zRange[0]  # Redshift
  a = 1./(1.+z)  # Scale factor
  rho_m = Mpart*(Ni/Lbox)**3  # Background density [(h/Mpc)^3]

  # Input files
  voidDir = voidBaseDir+'/'+sampleDir+'/'
  voidFilename1 = 'centers_central_'+sample.fullName+'.out'
  voidFilename2 = 'voidDesc_central_'+sample.fullName+'.out'


  # Read files
  matter_file = open(matterDir+matterFilename,'r')
  matter_data = []
  count = 0
  sep = 1e8
  for (i,line) in enumerate(matter_file):
    if i/int(sep) == i/sep: print str(round(i*100./Ni**3,1))+' percent of all particles read'
    if random.random() > ss: continue
  count += 1
  matter_data.append(line.split(',')[0:3])

  print str(count)+' particles read in total'
  matter_data = np.asarray(matter_data,dtype=np.float32)
  matter_file.close()


  halo_file = open(haloDir+haloFilename,'r')
  halo_data = np.reshape(halo_file.read().replace('\n',',').split(',')[0:-1],(-1,12)).astype(np.float32)
  halo_file.close()

  void_file = open(voidDir+voidFilename1,'r')
  void_header1 = void_file.readline().split(",")
  void_data1 = np.reshape(void_file.read().split(),(-1,len(void_header))).astype(np.float32)
  void_file.close()

  void_file = open(voidDir+voidFilename2,'r')
  void_header2 = void_file.readline().split(",")+void_file.readline().split(",")
  void_data2 = np.reshape(void_file.read().split(),(-1,11)).astype(np.float32)
  void_file.close()


  # Define arrays
  xm = matter_data[:,0:3]

  xh = halo_data[:,0:3]
  mh = halo_data[:,6]

  xv = void_data1[:,0:3]
  rv = void_data1[:,4]
  vv = void_data1[:,6]
  mv = void_data2[:,8]*Mpart


  # Interpolate to mesh
  dm, wm, ws = xcorlib.cic(xm, Lbox, Lboxcut = Lboxcut, Nmesh = Nmesh)
  dh, wm, ws = xcorlib.cic(xh, Lbox, Lboxcut = Lboxcut, Nmesh = Nmesh)
  dv, wm, ws = xcorlib.cic(xv, Lbox, Lboxcut = Lboxcut, Nmesh = Nmesh)

  # Load dark matter grid
  #output = open('dm_'+str(Nmesh)+'_ss'+str(ss)+'_z'+str(z)+'.dat', 'rb')
  #dm = pickle.load(output)
  #output.close()

  # Save dark matter grid
  #output = open('dm_'+str(Nmesh)+'_ss'+str(ss)+'_z'+str(z)+'.dat', 'wb')
  #pickle.dump(dm,output)
  #output.close()


  # Power spectra & correlation functions
  ((Nm, km, Pmm, SPmm),(Nmx, rm, Xmm, SXmm)) = xcorlib.powcor(dm, dm, Lbox, Nmesh = Nmesh, Nbin = Nbin, scale = 'lin', cor = True)
  ((Nm, km, Pvm, SPvm),(Nmx, rm, Xvm, SXvm)) = xcorlib.powcor(dv, dm, Lbox, Nmesh = Nmesh, Nbin = Nbin, scale = 'lin', cor = True)
  ((Nm, km, Phm, SPhm),(Nmx, rm, Xhm, SXhm)) = xcorlib.powcor(dh, dm, Lbox, Nmesh = Nmesh, Nbin = Nbin, scale = 'lin', cor = True)
  ((Nm, km, Pvv, SPvv),(Nmx, rm, Xvv, SXvv)) = xcorlib.powcor(dv, dv, Lbox, Nmesh = Nmesh, Nbin = Nbin, scale = 'lin', cor = True)
  ((Nm, km, Pvh, SPvh),(Nmx, rm, Xvh, SXvh)) = xcorlib.powcor(dv, dh, Lbox, Nmesh = Nmesh, Nbin = Nbin, scale = 'lin', cor = True)
  ((Nm, km, Phh, SPhh),(Nmx, rm, Xhh, SXhh)) = xcorlib.powcor(dh, dh, Lbox, Nmesh = Nmesh, Nbin = Nbin, scale = 'lin', cor = True)


  # Number densities
  nm = np.empty(len(km))
  nh = np.empty(len(km))
  nv = np.empty(len(km))
  nm[:] = Npart/Lbox**3
  nh[:] = len(xh)/Lbox**3
  nv[:] = len(xv)/Lbox**3

  # Number functions
  Mbin = 40
  Vbin = 40
  Nh, Mh = np.histogram(mh, bins = mh.min()*(mh.max()/mh.min())**(np.arange(Mbin+1)/float(Mbin)))
  Nvm, Mv = np.histogram(mv, bins = mv.min()*(mv.max()/mv.min())**(np.arange(Vbin+1)/float(Vbin)))
  Nvv, Vv = np.histogram(vv, bins = vv.min()*(vv.max()/vv.min())**(np.arange(Vbin+1)/float(Vbin)))

  # Bias
  b_hh = np.sqrt(Phh/Pmm)
  b_vv = np.sqrt(Pvv/Pmm)
  b_hm = Phm/Pmm
  b_vm = Pvm/Pmm
  b_vh = Pvh/Phh

  knl = 0.04  # Wavenumber above which nonlinearities kick in [h/Mpc]
  idls = np.where(km <= knl)[0]
  bm_hm = np.average(b_hm[idls],weights=Nm[idls])
  bm_vm = np.average(b_vm[idls],weights=Nm[idls])
  bm_vh = np.average(b_vh[idls],weights=Nm[idls])

  # Shot Noise
  sn_hh = Phh - Phm**2/Pmm
  sn_vh = Pvh - Pvm*Phm/Pmm
  sn_vv = Pvv - Pvm**2/Pmm



  # Plots
  ms = 4
  mew = 0.2
  margin = 1.2
  kmin = km.min()/margin
  kmax = km.max()*margin
  rmin = rm.min()/margin
  rmax = rm.max()*margin

  plt.imshow(np.sum(dm+1,2)/Nmesh,extent=[0,Lbox,0,Lbox],aspect='equal',cmap='YlGnBu_r',interpolation='gaussian')
  plt.xlabel(r'$x \;[h^{-1}\mathrm{Mpc}]$')
  plt.ylabel(r'$y \;[h^{-1}\mathrm{Mpc}]$')
  plt.title(r'Dark matter')
  plt.colorbar()
  plt.savefig(outputDir+'/dm_'+sample.fullName+'.pdf', format='pdf', bbox_inches="tight")
  plt.clf()

  plt.imshow(np.sum(dv+1,2)/Nmesh,extent=[0,Lbox,0,Lbox],aspect='equal',cmap='YlGnBu_r',interpolation='gaussian')
  plt.xlabel(r'$x \;[h^{-1}\mathrm{Mpc}]$')
  plt.ylabel(r'$y \;[h^{-1}\mathrm{Mpc}]$')
  plt.title(r'Voids')
  plt.colorbar()
  plt.savefig(outputDir+'/dv_'+sample.fullName+'.pdf', format='pdf', bbox_inches="tight")
  plt.clf()

  plt.imshow(np.sum(dh+1,2)/Nmesh,extent=[0,Lbox,0,Lbox],aspect='equal',cmap='YlGnBu_r',interpolation='gaussian')
  plt.xlabel(r'$x \;[h^{-1}\mathrm{Mpc}]$')
  plt.ylabel(r'$y \;[h^{-1}\mathrm{Mpc}]$')
  plt.title(r'Halos')
  plt.colorbar()
  plt.savefig(outputDir+'/dh_'+sample.fullName+'.pdf', format='pdf', bbox_inches="tight")
  plt.clf()


  pa, = plt.plot(Mh[:-1], Nh, 'ro-', ms=ms, mew=mew)
  pb, = plt.plot(Vv[:-1]*1e9, Nvv, 'bo-', ms=ms, mew=mew)
  plt.xlabel(r'$M \;[h^{-1}M_{\odot}]$ , $V \;[h^{-3}\mathrm{kpc}^3]$')
  plt.ylabel(r'$N(M,V)$')
  plt.title(r'Number of halos and voids')
  plt.xscale('log')
  plt.yscale('log')
  plt.xlim(min(10**np.floor(np.log10(Vv.min()))*1e9,10**np.floor(np.log10(Mh.min()))), max(10**np.ceil(np.log10(Mh.max())),10**np.ceil(np.log10(Vv.max()))*1e8))
  plt.ylim(10**np.floor(np.log10(Nh.min())), 10**np.ceil(np.log10(Nh.max())))
  plt.annotate(r'$\frac{4\pi}{3}\langle r_\mathrm{v}\rangle^3$', xy=(4./3.*np.pi*rv.mean()**3*1e9,1.1), xytext=(-50,235),textcoords='offset points',arrowprops=dict(fc='k',arrowstyle="->",connectionstyle="angle,angleA=0,angleB=90,rad=10"))
  plt.legend([pa,pb],['halos','voids'],'best' )
  plt.savefig(outputDir+'/number_'+sample.fullName+'.pdf', format='pdf', bbox_inches="tight")
  plt.clf()


  plt.subplot(211)
  plt.subplots_adjust(wspace=0,hspace=0)
  plt.plot(np.sort(vv), mv[np.argsort(vv)], 'r-', ms=ms, mew=mew, lw=0.01)
  plt.plot(np.sort(vv), np.sort(mv), 'k-', ms=ms, mew=mew)
  plt.plot(np.sort(vv), np.sort(vv)*1e9, 'k--', ms=ms, mew=mew)
  plt.title(r'Mass-volume relation of voids')
  plt.ylabel(r'$M \;[h^{-1}M_\odot]$')
  plt.xscale('log')
  plt.yscale('log')
  plt.subplot(211).xaxis.set_major_formatter(NullFormatter())
  plt.subplot(212)
  plt.plot(np.sort(vv), mv[np.argsort(vv)]/np.sort(vv)/rho_m-1., 'b-', ms=ms, mew=mew, lw=0.01)
  plt.plot(np.sort(vv), np.sort(mv)/np.sort(vv)/rho_m-1., 'k-', ms=ms, mew=mew)
  plt.xlabel(r'$V \;[h^{-3}\mathrm{Mpc}^3]$')
  plt.ylabel(r'$\delta$')
  plt.xscale('log')
  plt.yscale('linear')
  plt.ylim(-1.01,-0.861)
  plt.savefig(outputDir+'/massvol_'+sample.fullName+'.pdf', format='pdf', bbox_inches="tight")
  plt.clf()


  pa ,= plt.plot(km, Phh, 'r-', ms=ms, mew=mew)
  plt.plot(km, Phh-sn_hh, 'r:', ms=ms, mew=mew)
  plt.fill_between(km, Phh+SPhh, abs(Phh-SPhh), color='r', alpha=0.2)
  pb ,= plt.plot(km, Phm, 'y-', ms=ms, mew=mew)
  plt.fill_between(km, Phm+SPhm, abs(Phm-SPhm), color='y', alpha=0.2)
  pc ,= plt.plot(km, Pmm, 'k-', ms=ms, mew=mew)
  plt.plot(km, Pmm-1./nm, 'k:', ms=ms, mew=mew)
  plt.fill_between(km, Pmm+SPmm, abs(Pmm-SPmm), color='k', alpha=0.2)
  pd ,= plt.plot(km, Pvh, 'g-', ms=ms, mew=mew)
  plt.plot(km, -Pvh, 'g--', ms=ms, mew=mew)
  plt.plot(km, abs(Pvh-sn_vh), 'g:', ms=ms, mew=mew)
  plt.fill_between(km, abs(Pvh+SPvh), abs(Pvh-SPvh), color='g', alpha=0.2)
  pe ,= plt.plot(km, Pvm, 'm-', ms=ms, mew=mew)
  plt.plot(km, -Pvm, 'm--', ms=ms, mew=mew)
  plt.fill_between(km, abs(Pvm+SPvm), abs(Pvm-SPvm), color='m', alpha=0.2)
  pf ,= plt.plot(km, Pvv, 'b-', ms=ms, mew=mew)
  plt.plot(km, Pvv-sn_vv, 'b:', ms=ms, mew=mew)
  plt.fill_between(km, Pvv+SPvv, abs(Pvv-SPvv), color='b', alpha=0.2)
  plt.annotate(r'$\frac{\pi}{\langle r_\mathrm{v}\rangle}$', xy=(np.pi/(rv.mean()),1.01*10**np.floor(np.log10(abs(Pvh).min()))/margin), xytext=(10,280),textcoords='offset points',arrowprops=dict(fc='k',arrowstyle="->",connectionstyle="angle,angleA=0,angleB=90,rad=10"))
  plt.xlabel(r'$k \;[h\mathrm{Mpc}^{-1}]$')
  plt.ylabel(r'$P(k) \;[h^{-3}\mathrm{Mpc}^3]$')
  plt.title(r'Power spectra')
  plt.xscale('log')
  plt.yscale('log')
  plt.xlim(kmin,kmax)
  plt.ylim(10**np.floor(np.log10(abs(Pvh).min()))/margin, max(10**np.ceil(np.log10(Phh.max())),10**np.ceil(np.log10(Pvv.max())))*margin)
  plt.legend([pa, pb, pc, pd, pe, pf],['hh', 'hm', 'mm', 'vh', 'vm', 'vv'],'lower left' )
  plt.savefig(outputDir+'/power_'+sample.fullName+'.pdf', format='pdf', bbox_inches="tight") 
  plt.clf()

  
  pa ,= plt.plot(rm, Xhh, 'r-', ms=ms, mew=mew)
  plt.fill_between(rm, abs(Xhh+SXhh), abs(Xhh-SXhh), color='r', alpha=0.2)
  pb ,= plt.plot(rm, Xhm, 'y-', ms=ms, mew=mew)
  plt.fill_between(rm, abs(Xhm+SXhm), abs(Xhm-SXhm), color='y', alpha=0.2)
  pc ,= plt.plot(rm, Xmm, 'k-', ms=ms, mew=mew)
  plt.fill_between(rm, abs(Xmm+SXmm), abs(Xmm-SXmm), color='k', alpha=0.2)
  pd ,= plt.plot(rm, Xvh, 'g-', ms=ms, mew=mew)
  plt.plot(rm, -Xvh, 'g--', ms=ms, mew=mew)
  plt.fill_between(rm, abs(Xvh+SXvh), abs(Xvh-SXvh), color='g', alpha=0.2)
  pe ,= plt.plot(rm, Xvm, 'm-', ms=ms, mew=mew)
  plt.plot(rm, -Xvm, 'm--', ms=ms, mew=mew)
  plt.fill_between(rm, abs(Xvm+SXvm), abs(Xvm-SXvm), color='m', alpha=0.2)
  pf ,= plt.plot(rm, Xvv, 'b-', ms=ms, mew=mew)
  plt.fill_between(rm, abs(Xvv+SXvv), abs(Xvv-SXvv), color='b', alpha=0.2)
  plt.annotate(r'$\langle r_\mathrm{v}\rangle$', xy=(rv.mean(),1.01*10**np.floor(np.log10(abs(Xvh).min()))/margin), xytext=(10,300),textcoords='offset points',arrowprops=dict(fc='k',arrowstyle="->",connectionstyle="angle,angleA=0,angleB=90,rad=10"))
  plt.xlabel(r'$r \;[h^{-1}\mathrm{Mpc}]$')
  plt.ylabel(r'$\xi(r)$')
  plt.title(r'Correlation functions')
  plt.xscale('log')
  plt.yscale('log')
  plt.xlim(rmin,rmax)
  plt.ylim(10**np.floor(np.log10(abs(Xvh).min()))/margin, max(10**np.ceil(np.log10(Xhh.max())),10**np.ceil(np.log10(Xvv.max())))*margin)
  plt.legend([pa, pb, pc, pd, pe, pf],['hh', 'hm', 'mm', 'vh', 'vm', 'vv'],'lower left' )
  plt.savefig(outputDir+'/correlation_'+sample.fullName+'.pdf', format='pdf', bbox_inches="tight") 
  plt.clf()
  
  
  pa, = plt.plot(km, b_hh, 'r-', ms=ms, mew=mew)
  pb, = plt.plot(km, b_hm, 'r--', ms=ms, mew=mew)
  pc, = plt.plot(km, b_vv, 'b-', ms=ms, mew=mew)
  pd, = plt.plot(km, b_vm, 'b--', ms=ms, mew=mew)
  pe, = plt.plot(km, b_vh/bm_vh, 'g-', ms=ms, mew=mew)
  plt.plot(km, np.sin(km*rv.mean())/(km*rv.mean()), 'k:', ms=ms, mew=mew)
  plt.xlabel(r'$k \;[h\mathrm{Mpc}^{-1}]$')
  plt.ylabel(r'$b(k)$')
  plt.title(r'Bias')
  plt.xscale('log')
  plt.yscale('linear')
  plt.xlim(kmin,kmax)
  plt.ylim(np.floor(b_vm.min()),np.ceil(max(b_hh.max(),b_vv.max())))
  plt.legend([pa,pb,pc,pd,pe],['hh', 'hm', 'vv', 'vm', r'$\bar{u}_\mathrm{v}(k)$'],'best' )
  plt.savefig(outputDir+'/bias_'+sample.fullName+'.pdf', format='pdf', bbox_inches="tight")
  plt.clf()
  
  
  pa, = plt.plot(km, sn_hh, 'r-', ms=ms, mew=mew)
  pb, = plt.plot(km, sn_vh, 'g-', ms=ms, mew=mew)
  pc, = plt.plot(km, sn_vv, 'b-', ms=ms, mew=mew)
  plt.plot(km, abs(sn_vh), 'g:')
  pd, = plt.plot(km, 1/nh, 'r--')
  pe, = plt.plot(km, 1/nv, 'b-.')
  plt.xlabel(r'$k \;[h\mathrm{Mpc}^{-1}]$')
  plt.ylabel(r'$\sigma^2(k)$')
  plt.title(r'Shotnoise')
  plt.xscale('log')
  plt.yscale('log')
  plt.xlim(kmin,kmax)
  plt.ylim(10**np.floor(np.log10(abs(sn_vh).min())), 10**np.ceil(np.log10(sn_vv.max())))
  plt.legend([pa,pb,pc,pd,pe],['hh', 'vh', 'vv', r'$\bar{n}_\mathrm{h}^{-1}$', r'$\bar{n}_\mathrm{v}^{-1}$'],'best' )
  plt.savefig(outputDir+'/shotnoise_'+sample.fullName+'.pdf', format='pdf', bbox_inches="tight")
  plt.clf()
