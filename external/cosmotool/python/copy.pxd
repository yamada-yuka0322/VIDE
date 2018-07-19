cimport cython
cimport numpy as npx

ctypedef fused sum_type:
    cython.int
    cython.float
    npx.uint64_t
    npx.uint32_t

@cython.boundscheck(False)
cdef inline sum_type _mysum(sum_type[:] jobs) nogil:
  cdef sum_type s
  cdef npx.uint64_t N
  cdef int i

  s = 0
  N = jobs.shape[0]
  for i in xrange(N):
    s += jobs[i]
  return s


