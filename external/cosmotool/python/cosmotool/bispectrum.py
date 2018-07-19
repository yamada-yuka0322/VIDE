import numpy as np

try:
    import cffi
    import os

    _ffi = cffi.FFI()
    _ffi.cdef("""

    void CosmoTool_compute_bispectrum(
        double *delta_hat, size_t Nx, size_t Ny, size_t Nz,
        size_t *Ntriangles,
        double* B, double delta_k, size_t Nk ) ;
    void CosmoTool_compute_powerspectrum(
        double *delta_hat, size_t Nx, size_t Ny, size_t Nz,
        size_t *Ncounts,
        double* P, double delta_k, size_t Nk );

    """);

    _pathlib = os.path.dirname(os.path.abspath(__file__))
    _lib = _ffi.dlopen(os.path.join(_pathlib,"_cosmo_bispectrum.so"))
except Exception as e:
    print(repr(e))
    raise RuntimeError("Failed to initialize _cosmo_bispectrum module")


def bispectrum(delta, delta_k, Nk, fourier=True):
    """bispectrum(delta, fourier=True)
    
    Args:
        * delta: a 3d density field, can be Fourier modes if fourier set to True
        
        
    Return:
        * A 3d array of the binned bispectrum
"""

    if len(delta.shape) != 3:
        raise ValueError("Invalid shape for delta")
    try:
        delta_k = float(delta_k)
        Nk = int(Nk)
    except:
        raise ValueError()

    if not fourier:
        delta = np.fft.rfftn(delta)        
    N1,N2,N3 = delta.shape
    rN3 = (N3-1)*2
    delta_hat_buf = np.empty((N1*N2*N3*2),dtype=np.double)
    delta_hat_buf[::2] = delta.real.ravel()
    delta_hat_buf[1::2] = delta.imag.ravel()
    
    size_size = _ffi.sizeof("size_t")
    if size_size == 4:
        triangle_buf = np.zeros((Nk,Nk,Nk),dtype=np.int32)
    elif size_size == 8:
        triangle_buf = np.zeros((Nk,Nk,Nk),dtype=np.int64)
    else:
        raise RuntimeError("Internal error, do not know how to map size_t")

    B_buf = np.zeros((Nk*Nk*Nk*2), dtype=np.double)
    
    _lib.CosmoTool_compute_bispectrum( \
        _ffi.cast("double *", delta_hat_buf.ctypes.data), \
        N1, N2, rN3, \
        _ffi.cast("size_t *", triangle_buf.ctypes.data), \
        _ffi.cast("double *", B_buf.ctypes.data), \
        delta_k,  \
        Nk)
    B_buf = B_buf.reshape((Nk,Nk,Nk,2))
    return triangle_buf, B_buf[...,0]+1j*B_buf[...,1]

def powerspectrum(delta, delta_k, Nk, fourier=True):
    """powerspectrum(delta, fourier=True)
    
    Args:
        * delta: a 3d density field, can be Fourier modes if fourier set to True
        
        
    Return:
        * A 3d array of the binned bispectrum
"""

    if len(delta.shape) != 3:
        raise ValueError("Invalid shape for delta")
    try:
        delta_k = float(delta_k)
        Nk = int(Nk)
    except:
        raise ValueError()

    if not fourier:
        delta = np.fft.rfftn(delta)        
    N1,N2,N3 = delta.shape
    delta_hat_buf = np.empty((N1*N2*N3*2),dtype=np.double)
    delta_hat_buf[::2] = delta.real.ravel()
    delta_hat_buf[1::2] = delta.imag.ravel()
    
    size_size = _ffi.sizeof("size_t")
    if size_size == 4:
        count_buf = np.zeros((Nk,),dtype=np.int32)
    elif size_size == 8:
        count_buf = np.zeros((Nk,),dtype=np.int64)
    else:
        raise RuntimeError("Internal error, do not know how to map size_t")

    B_buf = np.zeros((Nk,), dtype=np.double)
    
    _lib.CosmoTool_compute_powerspectrum( \
        _ffi.cast("double *", delta_hat_buf.ctypes.data), \
        N1, N2, N3, \
        _ffi.cast("size_t *", count_buf.ctypes.data), \
        _ffi.cast("double *", B_buf.ctypes.data), \
        delta_k,  \
        Nk)
    return count_buf, B_buf[...]

        
if __name__=="__main__":
    delta=np.zeros((16,16,16))
    delta[0,0,0]=1
    delta[3,2,1]=1
    b = powerspectrum(delta, 1, 16, fourier=False)
    a = bispectrum(delta, 1, 16, fourier=False)
    print(a[0].max())
