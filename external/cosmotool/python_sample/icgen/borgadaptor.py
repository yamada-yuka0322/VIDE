import numpy as np

def fourier_analysis(borg_vol):
  L = (borg_vol.ranges[1]-borg_vol.ranges[0])
  N = borg_vol.density.shape[0]
    
  return np.fft.rfftn(borg_vol.density)*(L/N)**3, L, N

def borg_upgrade_sampling(dhat, supersample):
    N = dhat.shape[0]
    N2 = N * supersample
    dhat_new = np.zeros((N2, N2, N2/2+1), dtype=np.complex128)

    hN = N/2
    dhat_new[:hN, :hN, :hN+1] = dhat[:hN, :hN, :]
    dhat_new[:hN, (N2-hN):N2, :hN+1] = dhat[:hN, hN:, :]
    dhat_new[(N2-hN):N2, (N2-hN):N2, :hN+1] = dhat[hN:, hN:, :]
    dhat_new[(N2-hN):N2, :hN, :hN+1] = dhat[hN:, :hN, :]

    return dhat_new, N2

def half_pixel_shift(borg, doshift=False):

    dhat,L,N = fourier_analysis(borg)
    if not doshift:
      return dhat, L

    return bare_half_pixel_shift(dhat, L, N)

def bare_half_pixel_shift(dhat, L, N, doshift=False):

#    dhat_new,N2 = borg_upgrade_sampling(dhat, 2)
#    d = (np.fft.irfftn(dhat_new)*(N2/L)**3)[1::2,1::2,1::2]
#    del dhat_new
#    dhat = np.fft.rfftn(d)*(L/N)**3
#    return dhat, L

#    dhat2 = np.zeros((N,N,N),dtype=np.complex128)
#    dhat2[:,:,:N/2+1] = dhat
#    dhat2[N:0:-1, N:0:-1, N:N/2:-1] = np.conj(dhat[1:,1:,1:N/2])
#    dhat2[0, N:0:-1, N:N/2:-1] = np.conj(dhat[0, 1:, 1:N/2])
#    dhat2[N:0:-1, 0, N:N/2:-1] = np.conj(dhat[1:, 0, 1:N/2])
#    dhat2[0,0,N:N/2:-1] = np.conj(dhat[0, 0, 1:N/2])

    ik = np.fft.fftfreq(N,d=L/N)*2*np.pi
    phi = 0.5*L/N*(ik[:,None,None]+ik[None,:,None]+ik[None,None,:(N/2+1)])
#    phi %= 2*np.pi
    phase = np.cos(phi)+1j*np.sin(phi)
    dhat = dhat*phase    
    dhat[N/2,:,:] = 0
    dhat[:,N/2,:] = 0
    dhat[:,:,N/2] = 0
    
    return dhat, L

