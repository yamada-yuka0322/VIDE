from libcpp cimport bool
from libcpp cimport string as cppstring
from libcpp.vector cimport vector as cppvector
from cython.parallel cimport prange
import numpy as np
cimport numpy as np
from cpython cimport PyObject, Py_INCREF
cimport cython

np.import_array()

cdef extern from "sys/types.h":
   ctypedef np.int64_t int64_t

cdef extern from "loadSimu.hpp" namespace "CosmoTool":

  cdef cppclass SimuData:
    np.float_t BoxSize
    np.float_t time
    np.float_t Hubble

    np.float_t Omega_M
    np.float_t Omega_Lambda
    np.int64_t TotalNumPart
    np.int64_t NumPart
    int64_t *Id
    float *Pos[3]
    float *Vel[3]
    int *type
    float *Mass

    bool noAuto

  cdef const int NEED_GADGET_ID 
  cdef const int NEED_POSITION
  cdef const int NEED_VELOCITY
  cdef const int NEED_TYPE
  cdef const int NEED_MASS

cdef extern from "loadGadget.hpp" namespace "CosmoTool":

  SimuData *loadGadgetMulti(const char *fname, int id, int flags, int gformat) nogil except +
  void cxx_writeGadget "CosmoTool::writeGadget" (const char * s, SimuData *data) except +

cdef extern from "safe_gadget.hpp":
  SimuData *loadGadgetMulti_safe(cppstring.string s, int flags, int gformat) nogil
  SimuData **alloc_simudata(int num) nogil
  void del_simudata(SimuData **d) nogil

cdef extern from "cppHelper.hpp":
  int customCosmotoolHandler() nogil

cdef extern from "loadRamses.hpp" namespace "CosmoTool":
  SimuData *loadRamsesSimu(const char *basename, int id, int cpuid, bool dp, int flags) except +customCosmotoolHandler

class PySimulationBase(object):
  """
  This is the base class to representation Simulation in CosmoTool/python.
  """
  
  def getPositions(self):
    """
    getPositions(self)
    
    Returns:
        A list of three arrays holding the positions of the particles. 
        The i-th element is the i-th coordinate of each particle.
        It may be None if the positions were not requested.
    """
    raise NotImplementedError("getPositions is not implemented")
    
  def getVelocities(self):
    """
    getVelocities(self)
    
    Returns:
        A list of three arrays holding the velocities of the particles. 
        The i-th element is the i-th coordinate of each particle.
        It may be None if the velocities were not requested.
    """
    raise NotImplementedError("getVelocities is not implemented")
    
  def getIdentifiers(self):
    """
    getIdentifiers(self)
    
    Returns:
        Returns an integer array that hold the unique identifiers of
        each particle. 
        It may be None if the identifiers were not requested.
    """
    raise NotImplementedError("getIdentifiers is not implemented")

  def getTypes(self):
    """
    getTypes(self)
    
    Returns:
        Returns an integer array that hold the type of
        each particle. 
        It may be None if the types were not requested.
    """
    raise NotImplementedError("getTypes is not implemented")

  def getOmega_M(self):
    """
    getOmega_M(self)
    
    Returns:
        the mean matter density in the simulation, with respect
        to the critical density.
    """ 
    raise NotImplementedError("getOmega_M is not implemented")
  
  def getOmega_Lambda(self):
    """
    getOmega_Lambda(self)
    
    Returns:
        the mean dark energy density in the simulation, with respect
        to the critical density.
    """ 
    raise NotImplementedError("getOmega_Lambda is not implemented")

  def getTime(self):
    """
    getTime(self)
    
    Returns:
        the time the snapshot was taken in the simulation. It can
        have various units depending on the file format.
    """
    raise NotImplementedError("getTime is not implemented")

  def getHubble(self):
    """
    getHubble(self)
    
    Returns:
        the hubble constant in unit of 100 km/s/Mpc
    """
    raise NotImplementedError("getHubble is not implemented")

  def getBoxsize(self):
    """
    getBoxsize(self)
    
    Returns:
        the size of the simulation box. The length unit is not fixed,
        though it is customary to have it in Mpc/h if the loader has
        access to the unit normalization.
    """
    raise NotImplementedError("getBoxsize is not implemented")

  def getMasses(self):
    """
    getMasses(self)
    
    Returns:
        an array with the masses of each particles, in unspecified unit that
        depend on the loader.
    """
    raise NotImplementedError("getMasses is not implemented")

cdef class Simulation:
    """
    Simulation()
    
    Class that directly manages internal loaded data obtained from a loader
    """

    cdef list positions
    cdef list velocities
    cdef object identifiers
    cdef object types
    cdef object masses

    cdef SimuData *data

    property BoxSize:
        def __get__(Simulation self):
            return self.data.BoxSize
        
    property time:
        def __get__(Simulation self):
            return self.data.time

    property Hubble:
        def __get__(Simulation self):
            return self.data.Hubble

    property Omega_M:
        def __get__(Simulation self):
            return self.data.Omega_M
        
    property Omega_Lambda:
        def __get__(Simulation self):
            return self.data.Omega_Lambda
        
    property positions:
        def __get__(Simulation self):
            return self.positions
        
    property velocities:
        def __get__(Simulation self):
            return self.velocities

    property identifiers:
        def __get__(Simulation self):
            return self.identifiers

    property types:
        def __get__(Simulation self):
            return self.types

    property masses:
        def __get__(Simulation self):
            return self.masses
        
    property numParticles:
        def __get__(Simulation self):
            return self.data.NumPart


    property totalNumParticles:
        def __get__(Simulation self):
            return self.data.TotalNumPart

    def __cinit__(Simulation self):
        self.data = <SimuData *>0
  
    def __dealloc__(Simulation self):
        if self.data != <SimuData *>0:
          del self.data


class PySimulationAdaptor(PySimulationBase):
  """
  PySimulationAdaptor(PySimulationBase_)
  
  This class is an adaptor for an internal type to the loader. It defines
  all the methods of PySimulationBase.
  
  Attributes:
    simu: a Simulation_ object
  """
  def __init__(self,sim):
    self.simu = sim

  def getBoxsize(self):
    return self.simu.BoxSize
    
  def getPositions(self):
    return self.simu.positions

  def getTypes(self):
    return self.simu.types
    
  def getVelocities(self):
    return self.simu.velocities
    
  def getIdentifiers(self):
    return self.simu.identifiers

  def getTime(self):
    return self.simu.time
    
  def getHubble(self):
    return self.simu.Hubble
    
  def getOmega_M(self):
    return self.simu.Omega_M
    
  def getOmega_Lambda(self):
    return self.simu.Omega_Lambda

  def getMasses(self):
    return self.simu.masses

cdef class ArrayWrapper:
    cdef void* data_ptr
    cdef np.uint64_t size
    cdef int type_array
    
    cdef set_data(self, np.uint64_t size, int type_array, void* data_ptr):
        """ Set the data of the array
            This cannot be done in the constructor as it must recieve C-level
            arguments.

            Args:
                size (int): Length of the array.
                data_ptr (void*): Pointer to the data
        """
        self.data_ptr = data_ptr
        self.size = size
        self.type_array = type_array
     
    def __array__(self):
        """ Here we use the __array__ method, that is called when numpy
tries to get an array from the object."""
        cdef np.npy_intp shape[1]
        
        shape[0] = <np.npy_intp> self.size
        # Create a 1D array, of length 'size'
        ndarray = np.PyArray_SimpleNewFromData(1, shape, self.type_array, self.data_ptr)
        return ndarray
     
    def __dealloc__(self):
        """ Frees the array. This is called by Python when all the
references to the object are gone. """
        pass
 
cdef object wrap_array(void *p, np.uint64_t s, int typ):
    cdef np.ndarray ndarray
    cdef ArrayWrapper wrapper

    wrapper = ArrayWrapper()
    wrapper.set_data(s, typ, p)
    ndarray = np.array(wrapper, copy=False)
    ndarray.base = <PyObject*> wrapper
    Py_INCREF(wrapper)
    
    return ndarray


cdef object wrap_float_array(float *p, np.uint64_t s):
    return wrap_array(<void *>p, s, np.NPY_FLOAT32)

cdef object wrap_int64_array(int64_t* p, np.uint64_t s):
    return wrap_array(<void *>p, s, np.NPY_INT64)

cdef object wrap_int_array(int* p, np.uint64_t s):
    return wrap_array(<void *>p, s, np.NPY_INT)

cdef object wrap_simudata(SimuData *data, int flags):
  cdef Simulation simu

  simu = Simulation()
  simu.data = data
  if flags & NEED_POSITION:
     simu.positions = [wrap_float_array(data.Pos[i], data.NumPart) for i in xrange(3)]
  else:
     simu.positions = None

  if flags & NEED_VELOCITY:
     simu.velocities = [wrap_float_array(data.Vel[i], data.NumPart) for i in xrange(3)]
  else:
     simu.velocities = None

  if flags & NEED_GADGET_ID:
     simu.identifiers = wrap_int64_array(data.Id, data.NumPart)
  else:
     simu.identifiers = None

  if flags & NEED_TYPE:
     simu.types = wrap_int_array(data.type, data.NumPart)
  else:
     simu.types = None

  if flags & NEED_MASS:
     simu.masses = wrap_float_array(data.Mass, data.NumPart)
  else:
     simu.masses = None

  return simu

def loadGadget(str filename, int snapshot_id, int gadgetFormat = 1, bool loadPosition = True, bool loadVelocity = False, bool loadId = False, bool loadType = False, bool loadMass=False):
    """loadGadget(filename, snapshot_id, gadgetFormat = 1, loadPosition=True, loadVelocity=False, loadId=False, loadType=False)

    This function loads Gadget-1 snapshot format.
    
    If snapshot_id is negative then the snapshot is considered not to be part of
    a set of snapshots written by different cpu. Otherwise the filename is modified
    to reflect the indicated snapshot_id.

    Arguments:
        filename (str): input filename
        snapshot_id (int): identifier of the gadget file if it is a multi-file snapshot
    
    Keyword arguments:
        loadPosition (bool): whether to load positions
        loadVelocity (bool): whether to load velocities
        loadId (bool): whether to load unique identifiers
        loadType (bool): whether to set types to particles
        loadMass (bool): whether to set the mass of particles
    
    Returns:
        an PySimulationAdaptor instance. 
"""

    cdef int flags
    cdef SimuData *data
    cdef Simulation simu
    cdef const char *filename_bs

    flags = 0
    if loadPosition:
        flags |= NEED_POSITION
    if loadVelocity:
        flags |= NEED_VELOCITY
    if loadId:
        flags |= NEED_GADGET_ID
    if loadType:
        flags |= NEED_TYPE
    if loadMass:
        flags |= NEED_MASS

    filename_b = bytes(filename, 'utf-8')
    filename_bs = filename_b
    with nogil:
      data = loadGadgetMulti(filename_bs, snapshot_id, flags, gadgetFormat)
    if data == <SimuData*>0:
        return None

    return PySimulationAdaptor(wrap_simudata(data, flags))

def loadParallelGadget(object filename_list, int gadgetFormat = 1, bool loadPosition = True, bool loadVelocity = False, bool loadId = False, bool loadType = False, bool loadMass=False):
  """loadParallelGadget(filename list, gadgetFormat=1, loadPosition=True, loadVelocity=False, loadId=False, loadType=False)

  Arguments:
    filename (list): a list or tuple of filenames to load in parallel
  
  Keyword arguments:
    loadPosition (bool): indicate to load positions
    loadVelocity (bool): indicate to load velocities
    loadId (bool): indicate to load id
    loadType (bool): indicate to load particle types
    
  Returns:
    It loads a gadget-1 snapshot and return a cosmotool.PySimulationBase_ object.
  """
  cdef int flags, i, num_files
  cdef list out_arrays
  cdef SimuData ** data
  cdef SimuData * local_data
  cdef Simulation simu
  cdef cppvector[cppstring.string] filenames

  flags = 0
  if loadPosition:
     flags |= NEED_POSITION
  if loadVelocity:
     flags |= NEED_VELOCITY
  if loadId:
     flags |= NEED_GADGET_ID
  if loadType:
     flags |= NEED_TYPE
  if loadMass:
     flags |= NEED_MASS

  num_files = len(filename_list)
  filenames.resize(num_files)
  data = alloc_simudata(num_files)
  for i,l in enumerate(filename_list):
    filenames[i] = l.encode('utf-8')
  
  with nogil:
    for i in prange(num_files):
      local_data = loadGadgetMulti_safe(filenames[i], flags, gadgetFormat)
      data[i] = local_data
#      data[i] = loadGadgetMulti(filenames[i].c_str(), -1, flags)
  
  out_arrays = []
  for i in xrange(num_files):
    if data[i] == <SimuData*>0:        
      out_arrays.append(None)
    else:
      out_arrays.append(PySimulationAdaptor(wrap_simudata(data[i], flags)))

  del_simudata(data)

  return out_arrays

def writeGadget(str filename, object simulation):
  """writeGadget(filename, simulation)
  
  This function attempts to write the content of the simulation object into 
  a file named `filename` using a Gadget-1 file format.
  
  Arguments:
    filename (str): output filename
    simulation (PySimulationBase): a simulation object
  """
  cdef SimuData simdata
  cdef np.ndarray[np.float32_t, ndim=1] pos, vel
  cdef np.ndarray[np.int64_t, ndim=1] ids
  cdef np.int64_t NumPart
  cdef int j
  
  if not isinstance(simulation,PySimulationBase):
    raise TypeError("Second argument must be of type SimulationBase")
  
  NumPart = simulation.positions[0].size
  simdata.noAuto = True
  
  for j in xrange(3):
    pos = simulation.getPositions()[j]
    vel = simulation.getVelocities()[j]
    
    if pos.size != NumPart or vel.size != NumPart:
      raise ValueError("Invalid number of particles")
    
    simdata.Pos[j] = <float *>pos.data
    simdata.Vel[j] = <float *>vel.data
  
  ids = simulation.getIdentifiers()
  simdata.Id = <int64_t *>ids.data
  simdata.BoxSize = simulation.getBoxsize()
  simdata.time = simulation.getTime()
  simdata.Hubble = simulation.getHubble()
  simdata.Omega_M = simulation.getOmega_M()
  simdata.Omega_Lambda = simulation.getOmega_Lambda()
  simdata.TotalNumPart = NumPart
  simdata.NumPart = NumPart

  cxx_writeGadget(filename, &simdata)

def loadRamses(str basepath, int snapshot_id, int cpu_id, bool doublePrecision = False, bool loadPosition = True, bool loadVelocity = False, bool loadId = False, bool loadMass = False):
  """ loadRamses(basepath, snapshot_id, cpu_id, doublePrecision = False, loadPosition = True, loadVelocity = False)
  Loads the indicated snapshot based on the cpu id, snapshot id and basepath. It is important to specify the correct precision in doublePrecision otherwise the loading will fail. There is no way of auto-detecting properly the precision of the snapshot file.
  
  Args:
    basepath (str): the base directory of the snapshot
    snapshot_id (int): the snapshot id
    cpu_id (int): the cpu id of the file to load
  
  Keyword args:
    doublePrecision (bool): By default it is False, thus singlePrecision
    loadPosition (bool): Whether to load positions
    loadVelocity (bool): Whether to load velocities
    loadId (bool): Whether to load identifiers
    loadMass (bool): Whether to load mass value 
    
  Returns:
    An object derived from PySimulationBase_. 
"""
  cdef int flags
  cdef SimuData *data
  cdef Simulation simu

  flags = 0
  if loadPosition:
     flags |= NEED_POSITION
  if loadVelocity:
     flags |= NEED_VELOCITY
  if loadId:
     flags |= NEED_GADGET_ID
  if loadMass:
     flags |= NEED_MASS

  encpath = basepath.encode('utf-8') 
  try:
    data = loadRamsesSimu(encpath, snapshot_id, cpu_id, doublePrecision, flags)
    if data == <SimuData*>0:
      return None
  except RuntimeError as e:
    raise RuntimeError(str(e) + ' (check the float precision in snapshot)')

  return PySimulationAdaptor(wrap_simudata(data, flags))
  
