from cpython cimport bool
from cython cimport view
from cython.parallel import prange, parallel
from libc.math cimport sin, cos, abs, floor, round, sqrt
import numpy as np
cimport numpy as npx
cimport cython

__all__=["fast_interp"]


@cython.boundscheck(False)
@cython.cdivision(True)
def fast_interp(xmin0, dx0, A0, y0, out0, beyond_val=np.nan):
  cdef double rq, q
  cdef int iq
  cdef long i, Asize, ysize
  cdef npx.float64_t xmin, dx
  cdef npx.float64_t[:] out
  cdef npx.float64_t[:] A, y
  cdef npx.float64_t beyond

  beyond=beyond_val
  xmin = xmin0
  dx = dx0
  A = A0
  y = y0
  ysize = y.size
  out = out0
  Asize = A.size

  if out.size != ysize:
    raise ValueError("out and y must have the same size")

  with nogil:
    for i in prange(ysize):
      q = (y[i] - xmin) / dx
      iq = int(floor(q))
      rq = (q-iq)
      if iq+1 >= Asize or iq < 0:
        out[i] = beyond
      else:
        out[i] = rq * A[iq+1] + (1-rq)*A[iq]
