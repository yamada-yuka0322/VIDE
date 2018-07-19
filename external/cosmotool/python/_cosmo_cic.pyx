from libcpp cimport bool
from libcpp cimport string as cppstring
from libcpp cimport vector as cppvector
import numpy as np
cimport numpy as np
from cpython cimport PyObject, Py_INCREF
cimport cython

np.import_array()

cdef extern from "cic.hpp" namespace "CosmoTool":

  ctypedef float CICType
  ctypedef float Coordinates[3]
  
  cdef cppclass CICParticles:
    float mass
    Coordinates coords

  cdef cppclass CICFilter:
  
    CICFilter(np.uint32_t resolution, double L) nogil
    
    void resetMesh() nogil
    void putParticles(CICParticles* particles, np.uint32_t N) nogil
    void getDensityField(CICType *& field, np.uint32_t& res) nogil


@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
def leanCic(float[:,:] particles, float L, int Resolution):
  cdef cppvector.vector[CICParticles] *p
  cdef CICFilter *cic
  cdef np.uint64_t i
  cdef CICType *field
  cdef np.uint32_t dummyRes
  cdef np.ndarray[np.float64_t, ndim=3] out_field
  cdef np.ndarray[np.float64_t, ndim=1] out_field0
  cdef np.float64_t[:] out_field_buf
  cdef np.uint64_t j

  cic = new CICFilter(Resolution, L)
  print("Reset mesh")
  cic.resetMesh()
  
  if particles.shape[1] != 3:
    raise ValueError("Particles must be Nx3 array")
  
  print("Inserting particles")
#  p = new cppvector.vector[CICParticles](particles.shape[0])
#  for i in xrange(particles.shape[0]):
#    *p[i].mass = 1
#    *p[i].coords[0] = particles[i,0]
#    *p[i].coords[1] = particles[i,1]
#    *p[i].coords[2] = particles[i,2]

#  cic.putParticles(&p[0], particles.shape[0])
  del p
    
  print("Done")
  field = <CICType*>0
  dummyRes = 0
  cic.getDensityField(field, dummyRes)

  print("Got to allocate a numpy %dx%dx%d" % (dummyRes, dummyRes,dummyRes))
  
  out_field = np.empty((dummyRes, dummyRes, dummyRes), dtype=np.float64)
  out_field0 = out_field.reshape(out_field.size)
  out_field_buf = out_field
  print("Copy")
  for j in xrange(out_field_buf.size):
    out_field_buf[j] = field[j]
  
  del cic
  return out_field
