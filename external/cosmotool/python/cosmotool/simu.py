import warnings
from _cosmotool import *

class SimulationBare(PySimulationBase):

  def __init__(self, *args):
    if len(args) == 0:
      return
      
    if not isinstance(args[0], PySimulationBase):
      raise TypeError("Simulation object to mirror must be a PySimulationBase")
    s = args[0]

    self.positions = [q.copy() for q in s.getPositions()] if s.getPositions() is not None else None
    self.velocities = [q.copy() for q in s.getVelocities()] if s.getVelocities() is not None else None
    self.identifiers = s.getIdentifiers().copy() if s.getIdentifiers() is not None else None
    self.types = s.getTypes().copy() if s.getTypes() is not None else None
    self.boxsize = s.getBoxsize()
    self.time = s.getTime()
    self.Hubble = s.getHubble()
    self.Omega_M = s.getOmega_M()
    self.Omega_Lambda = s.getOmega_Lambda()
    try:
      self.masses = s.getMasses().copy() if s.getMasses() is not None else None
    except Exception as e:
      warnings.warn("Unexpected exception: " + repr(e))
  
  
  def merge(self, other):

    def _safe_merge(a, b):
      if b is not None:
        if a is not None:
          a = [np.append(q, r) for q,r in zip(a,b)]
        else:
          a = b
      return a
      
    def _safe_merge0(a, b):
      if b is not None:
        if a is not None:
          a = np.append(a, b)
        else:
          a = b
      return a
  
  
    assert self.time == other.getTime()
    assert self.Hubble == other.getHubble()
    assert self.boxsize == other.getBoxsize()
    assert self.Omega_M == other.getOmega_M()
    assert self.Omega_Lambda == other.getOmega_Lambda()
    
    self.positions = _safe_merge(self.positions, other.getPositions())
    self.velocities = _safe_merge(self.velocities, other.getVelocities())
    self.identifiers = _safe_merge0(self.identifiers, other.getIdentifiers())
    self.types = _safe_merge0(self.types, other.getTypes())
    try:
      self.masses = _safe_merge0(self.masses, other.getMasses())
    except Exception as e:
      warnings.warn("Unexpected exception: " + repr(e));
      self.masses = None

  def getTypes(self):
    return self.types
    
  def getPositions(self):
    return self.positions
    
  def getVelocities(self):
    return self.velocities
    
  def getIdentifiers(self):
    return self.identifiers

  def getMasses(self):
    return self.masses

  def getTime(self):
    return self.time
    
  def getHubble(self):
    return self.Hubble

  def getBoxsize(self):
    return self.boxsize

  def getOmega_M(self):
    return self.Omega_M

  def getOmega_Lambda(self):
    return self.Omega_Lambda


def simpleWriteGadget(filename, positions, boxsize=1.0, Hubble=100, Omega_M=0.30, time=1, velocities=None, identifiers=None):

  s = SimulationBare()
  
  s.positions = positions
  
  if velocities:
    s.velocities = velocities
  else:
    s.velocities = [np.zeros(positions[0].size,dtype=np.float32)]*3
    
  if identifiers:
    s.identifiers = identifiers
  else:
    s.identifiers = np.arange(positions[0].size, dtype=np.int64)
    
  s.Hubble = Hubble
  s.time = time
  s.Omega_M = Omega_M
  s.Omega_Lambda = 1-Omega_M
  s.boxsize = boxsize
    
  writeGadget(filename, s)
  
def loadRamsesAll(basepath, snapshot_id, **kwargs):
  """This function loads an entire ramses snapshot in memory. The keyword arguments are the one accepted
  by cosmotool.loadRamses
  
  Args:
    basepath (str): The base path of the snapshot (i.e. the directory holding the output_XXXXX directories)
    snapshot_id (int): The identifier of the snapshot to load.

  Keyword args:
    verboseMulti (bool): If true, print some progress information on loading multiple files
    
  See Also:
    cosmotool.loadRamses
    
  """
  cpu_id = 0
  output = None
  verbose = kwargs.get('verboseMulti',False)
  new_kw = dict(kwargs)
  if 'verboseMulti' in new_kw:
    del new_kw['verboseMulti']
  while True:
    base = "%s/output_%05d" % (basepath,snapshot_id)
    if verbose:
      print("Loading sub-snapshot %s (cpu_id=%d)" % (base,cpu_id))
    s = loadRamses(base, snapshot_id, cpu_id, **new_kw)
    if s == None:
      break
    if output == None:
      output = SimulationBare(s)
    else:
      output.merge(s)
      
    cpu_id += 1

  return output 
