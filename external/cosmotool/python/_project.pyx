from cpython cimport bool
from cython cimport view
from cython.parallel import prange, parallel
from libc.math cimport sin, cos, abs, floor, round, sqrt
import numpy as np
cimport numpy as npx
cimport cython
from copy cimport *

ctypedef npx.float64_t DTYPE_t
DTYPE=np.float64
FORMAT_DTYPE="d"

__all__=["project_cic","line_of_sight_projection","spherical_projection","DTYPE","interp3d","interp2d"]

cdef extern from "project_tool.hpp" namespace "":

  DTYPE_t compute_projection(DTYPE_t *vertex_value, DTYPE_t *u, DTYPE_t *u0, DTYPE_t rho) nogil

cdef extern from "openmp.hpp" namespace "CosmoTool":
  int smp_get_max_threads() nogil
  int smp_get_thread_id() nogil


@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
cdef void interp3d_INTERNAL_periodic(DTYPE_t x, DTYPE_t y, 
                                        DTYPE_t z, 
                                        DTYPE_t[:,:,:] d, DTYPE_t Lbox, DTYPE_t *retval) nogil:

    cdef int Ngrid = d.shape[0]
    cdef DTYPE_t inv_delta = Ngrid/Lbox
    cdef int ix, iy, iz
    cdef DTYPE_t f[2][2][2]
    cdef DTYPE_t rx, ry, rz
    cdef int jx, jy, jz

    rx = (inv_delta*x)
    ry = (inv_delta*y)
    rz = (inv_delta*z)

    ix = int(floor(rx))
    iy = int(floor(ry))
    iz = int(floor(rz))


    rx -= ix
    ry -= iy
    rz -= iz

    ix = ix % Ngrid
    iy = iy % Ngrid
    iz = iz % Ngrid

    jx = (ix+1)%Ngrid
    jy = (iy+1)%Ngrid
    jz = (iz+1)%Ngrid

    ix = ix%Ngrid
    iy = iy%Ngrid
    iz = iz%Ngrid

    f[0][0][0] = (1-rx)*(1-ry)*(1-rz)
    f[1][0][0] = (  rx)*(1-ry)*(1-rz)
    f[0][1][0] = (1-rx)*(  ry)*(1-rz)
    f[1][1][0] = (  rx)*(  ry)*(1-rz)

    f[0][0][1] = (1-rx)*(1-ry)*(  rz)
    f[1][0][1] = (  rx)*(1-ry)*(  rz)
    f[0][1][1] = (1-rx)*(  ry)*(  rz)
    f[1][1][1] = (  rx)*(  ry)*(  rz)

    retval[0] = \
        d[ix ,iy ,iz ] * f[0][0][0] + \
        d[jx ,iy ,iz ] * f[1][0][0] + \
        d[ix ,jy ,iz ] * f[0][1][0] + \
        d[jx ,jy ,iz ] * f[1][1][0] + \
        d[ix ,iy ,jz ] * f[0][0][1] + \
        d[jx ,iy ,jz ] * f[1][0][1] + \
        d[ix ,jy ,jz ] * f[0][1][1] + \
        d[jx ,jy ,jz ] * f[1][1][1]

@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
cdef void ngp3d_INTERNAL_periodic(DTYPE_t x, DTYPE_t y, 
                                 DTYPE_t z, 
                                 DTYPE_t[:,:,:] d, DTYPE_t Lbox, DTYPE_t *retval) nogil:

    cdef int Ngrid = d.shape[0]
    cdef DTYPE_t inv_delta = Ngrid/Lbox
    cdef int ix, iy, iz
    cdef DTYPE_t f[2][2][2]
    cdef DTYPE_t rx, ry, rz
    cdef int jx, jy, jz

    rx = (inv_delta*x)
    ry = (inv_delta*y)
    rz = (inv_delta*z)

    ix = int(round(rx))
    iy = int(round(ry))
    iz = int(round(rz))


    ix = ix%Ngrid
    iy = iy%Ngrid
    iz = iz%Ngrid

    retval[0] = d[ix ,iy ,iz ] 


@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
cdef void ngp3d_INTERNAL(DTYPE_t x, DTYPE_t y, 
                         DTYPE_t z, 
                         DTYPE_t[:,:,:] d, DTYPE_t Lbox, DTYPE_t *retval, DTYPE_t inval) nogil:

    cdef int Ngrid = d.shape[0]
    cdef DTYPE_t inv_delta = Ngrid/Lbox
    cdef int ix, iy, iz
    cdef DTYPE_t f[2][2][2]
    cdef DTYPE_t rx, ry, rz
    cdef int jx, jy, jz

    rx = (inv_delta*x)
    ry = (inv_delta*y)
    rz = (inv_delta*z)

    ix = int(round(rx))
    iy = int(round(ry))
    iz = int(round(rz))

    if ((ix < 0)  or (ix+1) >= Ngrid or (iy < 0)  or (iy+1) >= Ngrid or (iz < 0)  or (iz+1) >= Ngrid):
      retval[0] = inval
      return

    retval[0] = d[ix ,iy ,iz ] 


@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
cdef void interp3d_INTERNAL(DTYPE_t x, DTYPE_t y, 
                               DTYPE_t z, 
                               DTYPE_t[:,:,:] d, DTYPE_t Lbox, DTYPE_t *retval, DTYPE_t inval) nogil:
    
    cdef int Ngrid = d.shape[0]
    cdef DTYPE_t inv_delta = Ngrid/Lbox
    cdef int ix, iy, iz
    cdef DTYPE_t f[2][2][2]
    cdef DTYPE_t rx, ry, rz

    rx = (inv_delta*x)
    ry = (inv_delta*y)
    rz = (inv_delta*z)

    ix = int(floor(rx))
    iy = int(floor(ry))
    iz = int(floor(rz))

    rx -= ix
    ry -= iy
    rz -= iz

    if ((ix < 0)  or (ix+1) >= Ngrid or (iy < 0)  or (iy+1) >= Ngrid or (iz < 0)  or (iz+1) >= Ngrid):
      retval[0] = inval
      return
    #   assert ((ix >= 0) and ((ix+1) < Ngrid))
#    assert ((iy >= 0) and ((iy+1) < Ngrid))
#    assert ((iz >= 0) and ((iz+1) < Ngrid))

    f[0][0][0] = (1-rx)*(1-ry)*(1-rz)
    f[1][0][0] = (  rx)*(1-ry)*(1-rz)
    f[0][1][0] = (1-rx)*(  ry)*(1-rz)
    f[1][1][0] = (  rx)*(  ry)*(1-rz)

    f[0][0][1] = (1-rx)*(1-ry)*(  rz)
    f[1][0][1] = (  rx)*(1-ry)*(  rz)
    f[0][1][1] = (1-rx)*(  ry)*(  rz)
    f[1][1][1] = (  rx)*(  ry)*(  rz)

    retval[0] = \
        d[ix  ,iy  ,iz  ] * f[0][0][0] + \
        d[ix+1,iy  ,iz  ] * f[1][0][0] + \
        d[ix  ,iy+1,iz  ] * f[0][1][0] + \
        d[ix+1,iy+1,iz  ] * f[1][1][0] + \
        d[ix  ,iy  ,iz+1] * f[0][0][1] + \
        d[ix+1,iy  ,iz+1] * f[1][0][1] + \
        d[ix  ,iy+1,iz+1] * f[0][1][1] + \
        d[ix+1,iy+1,iz+1] * f[1][1][1]

@cython.boundscheck(False)
def interp3d(x not None, y not None, 
             z not None,
             npx.ndarray[DTYPE_t, ndim=3] d not None, DTYPE_t Lbox,
             bool periodic=False, bool centered=True, bool ngp=False, DTYPE_t inval = 0):
    """ interp3d(x,y,z,d,Lbox,periodic=False,centered=True,ngp=False) -> interpolated values
    
    Compute the tri-linear interpolation of the given field (d) at the given position (x,y,z). It assumes that they are box-centered coordinates. So (x,y,z) == (0,0,0) is equivalent to the pixel at (Nx/2,Ny/2,Nz/2) with Nx,Ny,Nz = d.shape. If periodic is set, it assumes the box is periodic 
"""
    cdef npx.ndarray[DTYPE_t] out
    cdef DTYPE_t[:] out_slice
    cdef DTYPE_t[:] ax, ay, az
    cdef DTYPE_t[:,:,:] in_slice
    cdef DTYPE_t retval
    cdef long i
    cdef long Nelt
    cdef int myperiodic, myngp
    cdef DTYPE_t shifter

    myperiodic = periodic
    myngp = ngp

    if centered:
      shifter = Lbox/2
    else:
      shifter = 0

    if d.shape[0] != d.shape[1] or d.shape[0] != d.shape[2]:
      raise ValueError("Grid must have a cubic shape")


    ierror = IndexError("Interpolating outside range")
    if type(x) == np.ndarray or type(y) == np.ndarray or type(z) == np.ndarray:
        if type(x) != np.ndarray or type(y) != np.ndarray or type(z) != np.ndarray:
            raise ValueError("All or no array. No partial arguments")
        
        ax = x
        ay = y
        az = z
        assert ax.size == ay.size and ax.size == az.size
        
        out = np.empty(x.shape, dtype=DTYPE)
        out_slice = out
        in_slice = d
        Nelt = ax.size
        with nogil:
          if not myngp:
            if myperiodic:
              for i in prange(Nelt):
                interp3d_INTERNAL_periodic(shifter+ax[i], shifter+ay[i], shifter+az[i], in_slice, Lbox, &out_slice[i])
            else:
              for i in prange(Nelt):
                interp3d_INTERNAL(shifter+ax[i], shifter+ay[i], shifter+az[i], in_slice, Lbox, &out_slice[i], inval)
          else:
            if myperiodic:
              for i in prange(Nelt):
                ngp3d_INTERNAL_periodic(shifter+ax[i], shifter+ay[i], shifter+az[i], in_slice, Lbox, &out_slice[i])
            else:
              for i in prange(Nelt):
                ngp3d_INTERNAL(shifter+ax[i], shifter+ay[i], shifter+az[i], in_slice, Lbox, &out_slice[i], inval)
        return out
    else:
      if not myngp:
        if periodic:
          interp3d_INTERNAL_periodic(shifter+x, shifter+y, shifter+z, d, Lbox, &retval)
        else:
          interp3d_INTERNAL(shifter+x, shifter+y, shifter+z, d, Lbox, &retval, inval)
      else:
        if periodic:
          ngp3d_INTERNAL_periodic(shifter+x, shifter+y, shifter+z, d, Lbox, &retval)
        else:
          ngp3d_INTERNAL(shifter+x, shifter+y, shifter+z, d, Lbox, &retval, inval)
      return retval

@cython.boundscheck(False)
@cython.cdivision(True)
cdef DTYPE_t interp2d_INTERNAL_periodic(DTYPE_t x, DTYPE_t y,
                                        npx.ndarray[DTYPE_t, ndim=2] d, DTYPE_t Lbox) except? 0:

    cdef int Ngrid = d.shape[0]
    cdef DTYPE_t inv_delta = Ngrid/Lbox
    cdef int ix, iy
    cdef DTYPE_t f[2][2]
    cdef DTYPE_t rx, ry
    cdef int jx, jy

    rx = (inv_delta*x + Ngrid/2)
    ry = (inv_delta*y + Ngrid/2)
    
    ix = int(floor(rx))
    iy = int(floor(ry))
    
    rx -= ix
    ry -= iy

    while ix < 0:
      ix += Ngrid
    while iy < 0:
      iy += Ngrid
    
    jx = (ix+1)%Ngrid
    jy = (iy+1)%Ngrid
    
    assert ((ix >= 0) and ((jx) < Ngrid))
    assert ((iy >= 0) and ((jy) < Ngrid))
    
    f[0][0] = (1-rx)*(1-ry)
    f[1][0] = (  rx)*(1-ry)
    f[0][1] = (1-rx)*(  ry)
    f[1][1] = (  rx)*(  ry)

    return \
        d[ix ,iy ] * f[0][0] + \
        d[jx ,iy ] * f[1][0] + \
        d[ix ,jy ] * f[0][1] + \
        d[jx ,jy ] * f[1][1]


@cython.boundscheck(False)
@cython.cdivision(True)
cdef DTYPE_t interp2d_INTERNAL(DTYPE_t x, DTYPE_t y,
                               npx.ndarray[DTYPE_t, ndim=2] d, DTYPE_t Lbox) except? 0:
    
    cdef int Ngrid = d.shape[0]
    cdef DTYPE_t inv_delta = Ngrid/Lbox
    cdef int ix, iy
    cdef DTYPE_t f[2][2]
    cdef DTYPE_t rx, ry

    rx = (inv_delta*x + Ngrid/2)
    ry = (inv_delta*y + Ngrid/2)

    ix = int(floor(rx))
    iy = int(floor(ry))

    rx -= ix
    ry -= iy

    if ((ix < 0)  or (ix+1) >= Ngrid):
      raise IndexError("X coord out of bound (ix=%d, x=%g)" % (ix,x))
    if ((iy < 0)  or (iy+1) >= Ngrid):
      raise IndexError("Y coord out of bound (iy=%d, y=%g)" % (iy,y))
#   assert ((ix >= 0) and ((ix+1) < Ngrid))
#    assert ((iy >= 0) and ((iy+1) < Ngrid))
#    assert ((iz >= 0) and ((iz+1) < Ngrid))

    f[0][0] = (1-rx)*(1-ry)
    f[1][0] = (  rx)*(1-ry)
    f[0][1] = (1-rx)*(  ry)
    f[1][1] = (  rx)*(  ry)

    return \
        d[ix  ,iy  ] * f[0][0] + \
        d[ix+1,iy  ] * f[1][0] + \
        d[ix  ,iy+1] * f[0][1] + \
        d[ix+1,iy+1] * f[1][1]
        
def interp2d(x not None, y not None,
             npx.ndarray[DTYPE_t, ndim=2] d not None, DTYPE_t Lbox,
             bool periodic=False):
    cdef npx.ndarray[DTYPE_t] out
    cdef npx.ndarray[DTYPE_t] ax, ay
    cdef int i

    if d.shape[0] != d.shape[1]:
      raise ValueError("Grid must have a square shape")

    if type(x) == np.ndarray or type(y) == np.ndarray:
        if type(x) != np.ndarray or type(y) != np.ndarray:
            raise ValueError("All or no array. No partial arguments")
        
        ax = x
        ay = y
        assert ax.size == ay.size 
        
        out = np.empty(x.shape, dtype=DTYPE)
        if periodic:
          for i in range(ax.size):
            out[i] = interp2d_INTERNAL_periodic(ax[i], ay[i], d, Lbox)
        else:
          for i in range(ax.size):
            out[i] = interp2d_INTERNAL(ax[i], ay[i], d, Lbox)

        return out
    else:
      if periodic:
        return interp2d_INTERNAL_periodic(x, y, d, Lbox)
      else:
        return interp2d_INTERNAL(x, y, d, Lbox)
        
 
@cython.boundscheck(False)
@cython.cdivision(True)
cdef void INTERNAL_project_cic_no_mass(DTYPE_t[:,:,:] g,
                                       DTYPE_t[:,:] x, int Ngrid, double Lbox, double shifter) nogil:
    cdef double delta_Box = Ngrid/Lbox
    cdef int i
    cdef double a[3]
    cdef double c[3]
    cdef int b[3]
    cdef int do_not_put

    for i in range(x.shape[0]):

        do_not_put = 0
        for j in range(3):
            a[j] = (x[i,j]+shifter)*delta_Box
            b[j] = int(floor(a[j]))
            a[j] -= b[j]
            c[j] = 1-a[j]
            if b[j] < 0 or b[j]+1 >= Ngrid:
              do_not_put = True

        if not do_not_put:
            g[b[0],b[1],b[2]] += c[0]*c[1]*c[2]
            g[b[0]+1,b[1],b[2]] += a[0]*c[1]*c[2]
            g[b[0],b[1]+1,b[2]] += c[0]*a[1]*c[2]
            g[b[0]+1,b[1]+1,b[2]] += a[0]*a[1]*c[2]

            g[b[0],b[1],b[2]+1] += c[0]*c[1]*a[2]
            g[b[0]+1,b[1],b[2]+1] += a[0]*c[1]*a[2]
            g[b[0],b[1]+1,b[2]+1] += c[0]*a[1]*a[2]
            g[b[0]+1,b[1]+1,b[2]+1] += a[0]*a[1]*a[2]

@cython.boundscheck(False)
@cython.cdivision(True)
cdef void INTERNAL_project_cic_no_mass_periodic(DTYPE_t[:,:,:] g,
                                                DTYPE_t[:,:] x, int Ngrid, double Lbox, double shifter) nogil:
    cdef double delta_Box = Ngrid/Lbox
    cdef int i
    cdef double a[3]
    cdef double c[3]
    cdef int b[3]
    cdef int b1[3]
    cdef int do_not_put
    cdef DTYPE_t[:,:] ax
    cdef DTYPE_t[:,:,:] ag

    ax = x
    ag = g

    for i in range(x.shape[0]):

        do_not_put = 0
        for j in range(3):
            a[j] = (ax[i,j]+shifter)*delta_Box
            b[j] = int(floor(a[j]))
            b1[j] = (b[j]+1) % Ngrid

            a[j] -= b[j]
            c[j] = 1-a[j]

            b[j] %= Ngrid

        ag[b[0],b[1],b[2]] += c[0]*c[1]*c[2]
        ag[b1[0],b[1],b[2]] += a[0]*c[1]*c[2]
        ag[b[0],b1[1],b[2]] += c[0]*a[1]*c[2]
        ag[b1[0],b1[1],b[2]] += a[0]*a[1]*c[2]
        
        ag[b[0],b[1],b1[2]] += c[0]*c[1]*a[2]
        ag[b1[0],b[1],b1[2]] += a[0]*c[1]*a[2]
        ag[b[0],b1[1],b1[2]] += c[0]*a[1]*a[2]
        ag[b1[0],b1[1],b1[2]] += a[0]*a[1]*a[2]


@cython.boundscheck(False)
@cython.cdivision(True)
cdef void INTERNAL_project_cic_with_mass(DTYPE_t[:,:,:] g,
                                         DTYPE_t[:,:] x,
                                         DTYPE_t[:] mass,
                                         int Ngrid, double Lbox, double shifter) nogil:
    cdef double delta_Box = Ngrid/Lbox
    cdef int i
    cdef double a[3]
    cdef double c[3]
    cdef DTYPE_t m0
    cdef int b[3]

    for i in range(x.shape[0]):

        do_not_put = False
        for j in range(3):
            a[j] = (x[i,j]+shifter)*delta_Box
            b[j] = int(a[j])
            a[j] -= b[j]
            c[j] = 1-a[j]
            if b[j] < 0 or b[j]+1 >= Ngrid:
              do_not_put = True

        if not do_not_put:
            m0 = mass[i]

            g[b[0],b[1],b[2]] += c[0]*c[1]*c[2]*m0
            g[b[0]+1,b[1],b[2]] += a[0]*c[1]*c[2]*m0
            g[b[0],b[1]+1,b[2]] += c[0]*a[1]*c[2]*m0
            g[b[0]+1,b[1]+1,b[2]] += a[0]*a[1]*c[2]*m0

            g[b[0],b[1],b[2]+1] += c[0]*c[1]*a[2]*m0
            g[b[0]+1,b[1],b[2]+1] += a[0]*c[1]*a[2]*m0
            g[b[0],b[1]+1,b[2]+1] += c[0]*a[1]*a[2]*m0
            g[b[0]+1,b[1]+1,b[2]+1] += a[0]*a[1]*a[2]*m0

@cython.boundscheck(False)
@cython.cdivision(True)
cdef void INTERNAL_project_cic_with_mass_periodic(DTYPE_t[:,:,:] g,
                                                  DTYPE_t[:,:] x,
                                                  DTYPE_t[:] mass,
                                                  int Ngrid, double Lbox, double shifter) nogil:
    cdef double half_Box = 0.5*Lbox, m0
    cdef double delta_Box = Ngrid/Lbox
    cdef int i
    cdef double a[3]
    cdef double c[3]
    cdef int b[3]
    cdef int b1[3]

    for i in range(x.shape[0]):

        for j in range(3):
            a[j] = (x[i,j]+shifter)*delta_Box
            b[j] = int(floor(a[j]))
            b1[j] = b[j]+1
            while b1[j] < 0:
              b1[j] += Ngrid
            while b1[j] >= Ngrid:
              b1[j] -= Ngrid

            a[j] -= b[j]
            c[j] = 1-a[j]

        m0 = mass[i]
        g[b[0],b[1],b[2]] += c[0]*c[1]*c[2]*m0
        g[b1[0],b[1],b[2]] += a[0]*c[1]*c[2]*m0
        g[b[0],b1[1],b[2]] += c[0]*a[1]*c[2]*m0
        g[b1[0],b1[1],b[2]] += a[0]*a[1]*c[2]*m0
        
        g[b[0],b[1],b1[2]] += c[0]*c[1]*a[2]*m0
        g[b1[0],b[1],b1[2]] += a[0]*c[1]*a[2]*m0
        g[b[0],b1[1],b1[2]] += c[0]*a[1]*a[2]*m0
        g[b1[0],b1[1],b1[2]] += a[0]*a[1]*a[2]*m0

       
def project_cic(npx.ndarray[DTYPE_t, ndim=2] x not None, npx.ndarray[DTYPE_t, ndim=1] mass, int Ngrid,
                double Lbox, bool periodic = False, centered=True):
    """
    project_cic(x array (N,3), mass (may be None), Ngrid, Lbox, periodict, centered=True)

    This function does a Cloud-In-Cell projection of a 3d unstructured dataset. First argument is a Nx3 array of coordinates. 
    Second argument is an optinal mass. Ngrid is the size output grid and Lbox is the physical size of the grid.
    """
    cdef npx.ndarray[DTYPE_t, ndim=3] g
    cdef double shifter
    cdef bool local_periodic

    local_periodic = periodic

    if centered:
      shifter = 0.5*Lbox
    else:
      shifter = 0

    if x.shape[1] != 3:
        raise ValueError("Invalid shape for x array")

    if mass is not None and mass.shape[0] != x.shape[0]:
        raise ValueError("Mass array and coordinate array must have the same number of elements")

    g = np.zeros((Ngrid,Ngrid,Ngrid),dtype=DTYPE)

    if not local_periodic:
        if mass is None:
         with nogil:
           INTERNAL_project_cic_no_mass(g, x, Ngrid, Lbox, shifter)
        else:
         with nogil:
           INTERNAL_project_cic_with_mass(g, x, mass, Ngrid, Lbox, shifter)       
    else:
        if mass is None:
         with nogil:
          INTERNAL_project_cic_no_mass_periodic(g, x, Ngrid, Lbox, shifter)
        else:
         with nogil:
          INTERNAL_project_cic_with_mass_periodic(g, x, mass, Ngrid, Lbox, shifter)       
         
    return g

def tophat_fourier_internal(npx.ndarray[DTYPE_t, ndim=1] x not None):
    cdef int i
    cdef npx.ndarray[DTYPE_t] y    
    cdef DTYPE_t x0

    y = np.empty(x.size, dtype=DTYPE)

    for i in range(x.size):
        x0 = x[i]
        if abs(x0)<1e-5:
            y[i] = 1
        else:
            y[i] = (3*(sin(x0) - x0 * cos(x0))/(x0**3))

    return y

def tophat_fourier(x not None):
    cdef npx.ndarray[DTYPE_t, ndim=1] b

    if type(x) != np.ndarray:
        raise ValueError("x must be a Numpy array")

    b = np.array(x, dtype=DTYPE).ravel()

    b = tophat_fourier_internal(b)

    return b.reshape(x.shape)

    

@cython.boundscheck(False)
@cython.cdivision(True)
cdef DTYPE_t cube_integral(DTYPE_t u[3], DTYPE_t u0[3], int r[1], DTYPE_t alpha_max) nogil:
    cdef DTYPE_t tmp_a
    cdef DTYPE_t v[3]
    cdef int i, j

    for i in xrange(3):
        if u[i] == 0.:
            continue

        if u[i] < 0:
            tmp_a = -u0[i]/u[i]
        else:
            tmp_a = (1-u0[i])/u[i]

        if tmp_a < alpha_max:
            alpha_max = tmp_a
            j = i

    for i in range(3):
        u0[i] += u[i]*alpha_max

    r[0] = j

    return alpha_max

@cython.boundscheck(False)
@cython.cdivision(True)
cdef DTYPE_t cube_integral_trilin(DTYPE_t u[3], DTYPE_t u0[3], int r[1], DTYPE_t vertex_value[8], DTYPE_t alpha_max) nogil:
    cdef DTYPE_t I, tmp_a
    cdef DTYPE_t v[3]
    cdef DTYPE_t term[4]
    cdef int i, j, q

    j = 0
    for i in range(3):
        if u[i] == 0.:
            continue

        if u[i] < 0:
            tmp_a = -u0[i]/u[i]
        else:
            tmp_a = (1-u0[i])/u[i]

        if tmp_a < alpha_max:
            alpha_max = tmp_a
            j = i
     
    I = compute_projection(vertex_value, u, u0, alpha_max)
    
    for i in xrange(3):
        u0[i] += u[i]*alpha_max

    # alpha_max is the integration length
    # we integrate between 0 and alpha_max (curvilinear coordinates)
    r[0] = j
    
    return I

@cython.boundscheck(False)
cdef DTYPE_t integrator0(DTYPE_t[:,:,:] density,
                         DTYPE_t u[3], DTYPE_t u0[3], int u_delta[3], int iu0[3], int jumper[1], DTYPE_t alpha_max) nogil:
    cdef DTYPE_t d
    
    d = density[iu0[0], iu0[1], iu0[2]]
    
    return cube_integral(u, u0, jumper, alpha_max)*d

@cython.boundscheck(False)
cdef DTYPE_t integrator1(DTYPE_t[:,:,:] density,
                         DTYPE_t u[3], DTYPE_t u0[3], int u_delta[3], int iu0[3], int jumper[1], DTYPE_t alpha_max) nogil:
    cdef DTYPE_t vertex_value[8]
    cdef DTYPE_t d
    cdef int a[3][2]
    cdef int i
    
    for i in xrange(3):
      a[i][0] = iu0[i]
      a[i][1] = iu0[i]+1

    vertex_value[0 + 2*0 + 4*0] = density[a[0][0], a[1][0], a[2][0]]
    vertex_value[1 + 2*0 + 4*0] = density[a[0][1], a[1][0], a[2][0]]
    vertex_value[0 + 2*1 + 4*0] = density[a[0][0], a[1][1], a[2][0]]
    vertex_value[1 + 2*1 + 4*0] = density[a[0][1], a[1][1], a[2][0]]

    vertex_value[0 + 2*0 + 4*1] = density[a[0][0], a[1][0], a[2][1]]
    vertex_value[1 + 2*0 + 4*1] = density[a[0][1], a[1][0], a[2][1]]
    vertex_value[0 + 2*1 + 4*1] = density[a[0][0], a[1][1], a[2][1]]
    vertex_value[1 + 2*1 + 4*1] = density[a[0][1], a[1][1], a[2][1]]

    return cube_integral_trilin(u, u0, jumper, vertex_value, alpha_max)


    
@cython.boundscheck(False)
cdef DTYPE_t C_line_of_sight_projection(DTYPE_t[:,:,:] density,
                             DTYPE_t a_u[3],
                             DTYPE_t min_distance,
                             DTYPE_t max_distance, DTYPE_t[:] shifter, int integrator_id) nogil except? 0:

    cdef DTYPE_t u[3] 
    cdef DTYPE_t ifu0[3]
    cdef DTYPE_t u0[3]
    cdef DTYPE_t utot[3]
    cdef int u_delta[3]
    cdef int iu0[3]
    cdef int i
    cdef int N = density.shape[0]
    cdef int half_N = density.shape[0]/2
    cdef int completed
    cdef DTYPE_t I0, d, dist2, delta, s, max_distance2
    cdef int jumper[1]
    
    cdef DTYPE_t (*integrator)(DTYPE_t[:,:,:],
                               DTYPE_t u[3], DTYPE_t u0[3], int u_delta[3], int iu0[3], int jumper[1], DTYPE_t alpha_max) nogil

    if integrator_id == 0:
      integrator = integrator0
    else:
      integrator = integrator1

    max_distance2 = max_distance**2

    for i in range(3):
        u[i] = a_u[i]
        u0[i] = a_u[i]*min_distance
        ifu0[i] = half_N+u0[i]+shifter[i]
        if (ifu0[i] <= 0 or ifu0[i] >= N):
          return 0
        iu0[i] = int(floor(ifu0[i]))
        u0[i] = ifu0[i]-iu0[i]
        u_delta[i] = 1 if iu0[i] > 0 else -1
        if (not ((iu0[i]>= 0) and (iu0[i] < N))):
          with gil:
            raise RuntimeError("iu0[%d] = %d !!" % (i,iu0[i]))
            
        if (not (u0[i]>=0 and u0[i]<=1)):
          with gil:
            raise RuntimeError("u0[%d] = %g !" % (i,u0[i]))

    completed = 0
    if ((iu0[0] >= N-1) or (iu0[0] <= 0) or
        (iu0[1] >= N-1) or (iu0[1] <= 0) or
        (iu0[2] >= N-1) or (iu0[2] <= 0)):
      completed = 1 

    I0 = 0
    jumper[0] = 0
    dist2 = 0
    while completed == 0:
        I0 += integrator(density, u, u0, u_delta, iu0, jumper, max_distance-sqrt(dist2))

        if u[jumper[0]] < 0:
            iu0[jumper[0]] -= 1
            u0[jumper[0]] = 1
        else:
            iu0[jumper[0]] += 1
            u0[jumper[0]] = 0

            
        if ((iu0[0] >= N-1) or (iu0[0] <= 0) or 
            (iu0[1] >= N-1) or (iu0[1] <= 0) or
            (iu0[2] >= N-1) or (iu0[2] <= 0)):
            completed = 1
        else:
          dist2 = 0
          for i in range(3):
             delta =  iu0[i]+u0[i]-half_N-shifter[i]
             dist2 += delta*delta

          if (dist2 > max_distance2):
            # Remove the last portion of the integral
            #delta = sqrt(dist2) - max_distance
            #I0 -= d*delta
            completed = 1
            
    return I0

def line_of_sight_projection(DTYPE_t[:,:,:] density not None,
                             DTYPE_t[:] a_u not None,
                             DTYPE_t min_distance,
                             DTYPE_t max_distance, DTYPE_t[:] shifter not None, int integrator_id=0):
  cdef DTYPE_t u[3]
  
  u[0] = a_u[0]
  u[1] = a_u[1]
  u[2] = a_u[2]
  
  return C_line_of_sight_projection(density,
                             u,
                             min_distance,
                             max_distance, shifter, integrator_id)

cdef double _spherical_projloop(double theta, double phi, DTYPE_t[:,:,:] density, 
                double min_distance, double max_distance, 
                DTYPE_t[:] shifter, int integrator_id) nogil:
  cdef DTYPE_t u0[3]

  stheta = sin(theta)
  u0[0] = cos(phi)*stheta
  u0[1] = sin(phi)*stheta
  u0[2] = cos(theta)
 
  return C_line_of_sight_projection(density, u0, min_distance, max_distance, shifter, integrator_id)


@cython.boundscheck(False)
def spherical_projection(int Nside,
                         npx.ndarray[DTYPE_t, ndim=3] density not None,
                         DTYPE_t min_distance,
                         DTYPE_t max_distance, int progress=1, int integrator_id=0, DTYPE_t[:] shifter = None, int booster=-1):
    """
    spherical_projection(Nside, density, min_distance, max_distance, progress=1, integrator_id=0, shifter=None, booster=-1)
    
    Keyword arguments:
      progress (int): show progress if it is equal to 1
      integrator_id (int): specify the order of integration along the line of shift
      shifter (DTYPE_t array): this is an array of size 3. It specifies the amount of shift to apply to the center, in unit of voxel
      booster (int): what is the frequency of refreshment of the progress bar. Small number decreases performance by locking the GIL.
      

    Arguments:
      Nside (int): Nside of the returned map
      density (NxNxN array): this is the density field, expressed as a cubic array
      min_distance (float): lower bound of the integration
      max_distance (float): upper bound of the integration
      
    Returns:
      an healpix map, as a 1-dimensional array.
    """
    import healpy as hp
    import progressbar as pb
    cdef int i
    cdef DTYPE_t[:] theta,phi
    cdef DTYPE_t[:,:,:] density_view
    cdef DTYPE_t[:] outm
    cdef int[:] job_done
    cdef npx.ndarray[DTYPE_t, ndim=1] outm_array
    cdef long N, N0
    cdef double stheta
    cdef int tid
    
    if shifter is None:
        shifter = view.array(shape=(3,), format=FORMAT_DTYPE, itemsize=sizeof(DTYPE_t))
        shifter[:] = 0
        
    print("allocating map")
    outm_array = np.empty(hp.nside2npix(Nside),dtype=DTYPE)
    print("initializing views")
    outm = outm_array
    density_view = density

    print("progress?")
    if progress != 0:
      p = pb.ProgressBar(maxval=outm.size,widgets=[pb.Bar(), pb.ETA()]).start()

    N = smp_get_max_threads()
    N0 = outm.size
    
    if booster < 0:
      booster = 1#000
    
    job_done = view.array(shape=(N,), format="i", itemsize=sizeof(int))
    job_done[:] = 0
    theta,phi = hp.pix2ang(Nside, np.arange(N0))
    with nogil, parallel():
      tid = smp_get_thread_id()
      for i in prange(N0,schedule='dynamic',chunksize=256):
          if progress != 0 and (i%booster) == 0:
            with gil:
              p.update(_mysum(job_done))
          outm[i] = _spherical_projloop(theta[i], phi[i], density_view, min_distance, max_distance, shifter, integrator_id)
          job_done[tid] += 1

    if progress:
      p.finish()


    return outm_array
