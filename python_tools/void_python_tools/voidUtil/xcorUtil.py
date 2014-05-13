import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
import xcorlib

def computeXcor(catalog,
                figDir="./",
                Nmesh = 256,
                Nbin = 100
                ):

# Computes and plots void-void and void-matter(galaxy) correlations
#   catalog: catalog to analyze
#   figDir: where to place plots
#   Nmesh: number of grid cells in power spectrum calculation
#   Nbin: number of grid cells in plot

  # Parameters
  Lbox = catalog.boxLen # Boxlength
  Lboxcut = 0.
  Lbox -= 2*Lboxcut
  
  # Input particle arrays of shape (N,3)
  xm = catalog.partPos # Halos / Galaxies / Dark matter
  xv = getArray(catalog.voids, 'barycenter')
  
  
  # Interpolate to mesh
  dm, wm, ws = xcorlib.cic(xm, Lbox, Lboxcut = Lboxcut, Nmesh = Nmesh, weights = None)
  dv, wm, ws = xcorlib.cic(xv, Lbox, Lboxcut = Lboxcut, Nmesh = Nmesh, weights = None)
  
  # Fourier transform
  dmk = np.fft.rfftn(dm)
  dvk = np.fft.rfftn(dv)
  
  # 1D Power spectra & correlation functions
  ((Nm, km, Pmm, SPmm),(Nmx, rm, Xmm, SXmm)) = xcorlib.powcor(dmk, dmk, Lbox, Nbin, 'lin', True, True, 1)
  ((Nm, km, Pvm, SPvm),(Nmx, rm, Xvm, SXvm)) = xcorlib.powcor(dvk, dmk, Lbox, Nbin, 'lin', True, True, 1)
  ((Nm, km, Pvv, SPvv),(Nmx, rm, Xvv, SXvv)) = xcorlib.powcor(dvk, dvk, Lbox, Nbin, 'lin', True, True, 1)
  
  # 2D Power spectra & correlation functions
  ((Nm2d, kmper, kmpar, Pmm2d),(Nmx2d, rmper, rmpar, Xmm2d)) = xcorlib.powcor(dmk, dmk, Lbox, Nbin, 'lin', True, True, 2)
  ((Nm2d, kmper, kmpar, Pvm2d),(Nmx2d, rmper, rmpar, Xvm2d)) = xcorlib.powcor(dvk, dmk, Lbox, Nbin, 'lin', True, True, 2)
  
  # Number densities
  nm = np.empty(len(km))
  nh = np.empty(len(km))
  nv = np.empty(len(km))
  nm[:] = len(xm)/Lbox**3
  nh[:] = len(xh)/Lbox**3
  nv[:] = len(xv)/Lbox**3
  
  # Bias
  b_hh = np.sqrt(Phh/Pmm)
  b_vv = np.sqrt(Pvv/Pmm)
  b_hm = Phm/Pmm
  b_vm = Pvm/Pmm
  b_vh = Pvh/Phh
  
  # Shot Noise
  sn_hh = Phh - Phm**2/Pmm
  sn_vh = Pvh - Pvm*Phm/Pmm
  sn_vv = Pvv - Pvm**2/Pmm
  
  
  # Plots
  mpl.rc('font', family='serif')
  ms = 2.5
  fs = 16
  mew = 0.1
  margin = 1.2
  kmin = km.min()/margin
  kmax = km.max()*margin
  rmin = rm.min()/margin
  rmax = rm.max()*margin
  
  
  # Density fields (projected)
  plt.imshow(np.sum(dm[:,:,:]+1,2),extent=[0,Lbox,0,Lbox],aspect='equal',cmap='YlGnBu_r',interpolation='gaussian')
  plt.xlabel(r'$x \;[h^{-1}\mathrm{Mpc}]$')
  plt.ylabel(r'$y \;[h^{-1}\mathrm{Mpc}]$')
  plt.title(r'Dark matter')
  plt.savefig(figDir+'/dm_'+sample.fullName+'.pdf', format='pdf', bbox_inches="tight")
  plt.clf()
  
  plt.imshow(np.sum(dv[:,:,:]+1,2)/Nmesh,extent=[0,Lbox,0,Lbox],aspect='equal',cmap='YlGnBu_r',interpolation='gaussian')
  plt.xlabel(r'$x \;[h^{-1}\mathrm{Mpc}]$')
  plt.ylabel(r'$y \;[h^{-1}\mathrm{Mpc}]$')
  plt.title(r'Voids')
  plt.savefig(figDir+'/dv_'+sample.fullName+'.pdf', format='pdf', bbox_inches="tight") #, dpi=300
  plt.clf()
  
  
  # Power spectra & correlation functions
  pa ,= plt.plot(km, Phh, 'r-s', ms=ms, mew=mew, mec='k')
  plt.plot(km, Phh-sn_hh, 'r--', ms=ms, mew=mew)
  plt.plot(km, sn_hh, 'r:', ms=ms, mew=mew)
  plt.fill_between(km, Phh+SPhh, abs(Phh-SPhh), color='r', alpha=0.2)
  pb ,= plt.plot(km, Phm, 'y-^', ms=ms, mew=mew, mec='k')
  plt.fill_between(km, Phm+SPhm, abs(Phm-SPhm), color='y', alpha=0.2)
  pc ,= plt.plot(km, Pmm, 'k-o', ms=0.8*ms, mew=mew, mec='k')
  plt.plot(km, Pmm-1./nm, 'k--', ms=ms, mew=mew)
  plt.fill_between(km, Pmm+SPmm, abs(Pmm-SPmm), color='k', alpha=0.2)
  pd ,= plt.plot(km, Pvh, 'g-*', ms=1.5*ms, mew=mew, mec='k')
  plt.plot(km, -Pvh, 'g*', ms=1.5*ms, mew=mew, mec='k')
  plt.plot(km, abs(Pvh-sn_vh), 'g--', ms=ms, mew=mew)
  plt.plot(km, sn_vh, 'g:', ms=ms, mew=mew)
  plt.plot(km, -sn_vh, 'g-.', ms=ms, mew=mew)
  plt.fill_between(km, abs(Pvh+SPvh), abs(Pvh-SPvh), color='g', alpha=0.2)
  pe ,= plt.plot(km, Pvm, 'm-D', ms=ms, mew=mew, mec='k')
  plt.plot(km, -Pvm, 'mD', ms=ms, mew=mew, mec='k')
  plt.fill_between(km, abs(Pvm+SPvm), abs(Pvm-SPvm), color='m', alpha=0.2)
  pf ,= plt.plot(km, Pvv, 'b-p', ms=1.3*ms, mew=mew, mec='k')
  plt.plot(km, Pvv-sn_vv, 'b--', ms=ms, mew=mew)
  plt.plot(km, sn_vv, 'b:', ms=ms, mew=mew)
  plt.fill_between(km, Pvv+SPvv, abs(Pvv-SPvv), color='b', alpha=0.2)
  plt.xlabel(r'$k \;[h\mathrm{Mpc}^{-1}]$')
  plt.ylabel(r'$P(k) \;[h^{-3}\mathrm{Mpc}^3]$')
  plt.title(r'Power spectra')
  plt.xscale('log')
  plt.yscale('log')
  plt.xlim(kmin,kmax)
  plt.ylim(10**np.floor(np.log10(abs(Pvh).min()))/margin, max(10**np.ceil(np.log10(Phh.max())),10**np.ceil(np.log10(Pvv.max())))*margin)
  plt.legend([pa, pb, pc, pd, pe, pf],['gg', 'gm', 'mm', 'vg', 'vm', 'vv'],'lower left',prop={'size':12})
  plt.savefig(figDir+'/power_'+sample.fullName+'.pdf', format='pdf', bbox_inches="tight") 
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
  plt.plot(rm, -Xvv, 'b--', ms=ms, mew=mew)
  plt.fill_between(rm, abs(Xvv+SXvv), abs(Xvv-SXvv), color='b', alpha=0.2)
  plt.xlabel(r'$r \;[h^{-1}\mathrm{Mpc}]$')
  plt.ylabel(r'$\xi(r)$')
  plt.title(r'Correlation functions')
  plt.xscale('log')
  plt.yscale('log')
  plt.xlim(rmin,rmax)
  plt.ylim(10**np.floor(np.log10(abs(Xvh).min()))/margin, max(10**np.ceil(np.log10(Xhh.max())),10**np.ceil(np.log10(Xvv.max())))*margin)
  plt.legend([pa, pb, pc, pd, pe, pf],['gg', 'gm', 'mm', 'vg', 'vm', 'vv'],'best',prop={'size':12})
  plt.savefig(figDir+'/correlation_'+sample.fullName+'.pdf', format='pdf', bbox_inches="tight") 
  plt.clf()
  
  
  # 2D power spectra & correlation functions
  kpermin = kmper.min()
  kpermax = 0.3001
  kparmin = kmpar.min()
  kparmax = 0.3001
  rpermin = rmper.min()
  rpermax = 40
  rparmin = rmpar.min()
  rparmax = 40
  
  for (P2d,idx,vmin,vmax) in ([Pmm2d,'mm',None,None],[Pvm2d,'vm',None,None],[Phm2d,'gm',None,None],[Pvv2d,'vv',None,2.9],[Pvh2d,'vg',None,None],[Phh2d,'gg',None,None]):
    cut = np.where(kmper[:,0] <= kpermax)[0].max()+2
    plt.pcolormesh(kmper[0:cut,0:cut], kmpar[0:cut,0:cut], P2d[0:cut,0:cut]/1e4, cmap=cm.Spectral_r, shading='gouraud', vmin=vmin, vmax=vmax)
    plt.colorbar(format='%.1f')
    plt.contour(kmper[0:cut,0:cut], kmpar[0:cut,0:cut], P2d[0:cut,0:cut]/1e4, levels=np.array(P2d.min()+(P2d.max()-P2d.min())*(np.arange(Nbin/6+1)/float(Nbin/6)))/1e4, vmin=vmin, vmax=vmax, colors='k', linewidths=0.2)
    plt.xlabel(r'$k_\perp \;[h\mathrm{Mpc}^{-1}]$', fontsize=fs)
    plt.ylabel(r'$k_\parallel \;[h\mathrm{Mpc}^{-1}]$', fontsize=fs)
    plt.axes().set_aspect('equal')
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlim(kpermin,kpermax)
    plt.ylim(kparmin,kparmax)
    plt.title(r'$P_{\mathrm{'+idx+r'}}(k_\perp, k_\parallel) \;[10^4h^{-3}\mathrm{Mpc}^3]$', fontsize=fs)
    plt.savefig(figDir+'/P'+idx+'2d_'+sample.fullName+'.pdf', format='pdf', bbox_inches="tight")
    plt.clf()
  
  for (X2d,idx,vmin,vmax) in ([Xmm2d,'mm',None,None],[Xvm2d,'vm',None,None],[Xhm2d,'gm',None,None],[Xvv2d,'vv',None,0.2],[Xvh2d,'vg',None,None],[Xhh2d,'gg',None,None]):
    cut = np.where(rmper[:,0] <= rpermax)[0].max()+3
    plt.pcolormesh(rmper[0:cut,0:cut], rmpar[0:cut,0:cut], X2d[0:cut,0:cut], cmap=cm.Spectral_r, shading='gouraud', vmin=vmin, vmax=vmax)
    plt.colorbar(format='%+.2f')
    plt.contour(rmper[0:cut,0:cut], rmpar[0:cut,0:cut], X2d[0:cut,0:cut], levels=np.array(X2d.min()+(X2d.max()-X2d.min())*(np.arange(Nbin/6+1)/float(Nbin/6))), vmin=vmin, vmax=vmax, colors='k', linewidths=0.2)  
    plt.xlabel(r'$r_\perp \;[h^{-1}\mathrm{Mpc}]$', fontsize=fs)
    plt.ylabel(r'$r_\parallel \;[h^{-1}\mathrm{Mpc}]$', fontsize=fs)
    plt.axes().set_aspect('equal')
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlim(rpermin,rpermax)
    plt.ylim(rparmin,rparmax)
    plt.title(r'$\xi_{\mathrm{'+idx+r'}}(r_\perp, r_\parallel)$', fontsize=fs)
    plt.savefig(figDir+'/X'+idx+'2d_'+sample.fullName+'.pdf', format='pdf', bbox_inches="tight")
    plt.clf()
 
  return 
