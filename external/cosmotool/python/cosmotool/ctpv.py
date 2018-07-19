import numpy as np
from contextlib import contextmanager

class ProgrammableParticleLoad(object):

  @staticmethod
  def main_script(source, particles, aname="default", aux=None):
    import vtk
    from vtk.util import numpy_support as ns
  
    out = source.GetOutput()
    vv = vtk.vtkPoints()

    assert len(particles.shape) == 2
    assert particles.shape[1] == 3

    vv.SetData(ns.numpy_to_vtk(np.ascontiguousarray(particles.astype(np.float64)), deep=1))
    vv.SetDataTypeToDouble()

    out.Allocate(1,1)
    out.SetPoints(vv)

    if aux is not None:
      for n,a in aux:
        a_vtk = ns.numpy_to_vtk(
           np.ascontiguousarray(a.astype(np.float64)
             ), 
           deep=1)
        a_vtk.SetName(n)
        out.GetPointData().AddArray(a_vtk)

    out.InsertNextCell(vtk.VTK_VERTEX, particles.shape[0], range(particles.shape[0]))

  @staticmethod
  def request_information(source):
   pass



class ProgrammableParticleHistoryLoad(object):

  @staticmethod
  def main_script(source, particles, velocities=None, aname="default",addtime=False):
    import vtk
    from vtk.util import numpy_support as ns
  
    out = source.GetOutput()
    vv = vtk.vtkPoints()

    assert len(particles.shape) == 3
    assert particles.shape[2] == 3
    
    if not velocities is None:
      for i,j in zip(velocities.shape,particles.shape):
        assert i==j

    Ntime,Npart,_ = particles.shape

    vv.SetData(ns.numpy_to_vtk(np.ascontiguousarray(particles.reshape((Ntime*Npart,3)).astype(np.float64)), deep=1))
    vv.SetDataTypeToDouble()

    out.Allocate(1,1)
    out.SetPoints(vv)

    if not velocities is None:
      print("Adding velocities")
      vel_vtk = ns.numpy_to_vtk(np.ascontiguousarray(velocities.reshape((Ntime*Npart,3)).astype(np.float64)), deep=1)
      vel_vtk.SetName("velocities")
      out.GetPointData().AddArray(vel_vtk)

    if addtime:
      timearray = np.arange(Ntime)[:,None].repeat(Npart, axis=1).reshape(Ntime*Npart)
      timearray = ns.numpy_to_vtk(np.ascontiguousarray(timearray.astype(np.float64)), deep=1)
      timearray.SetName("timearray")
      out.GetPointData().AddArray(timearray)

    out.InsertNextCell(vtk.VTK_VERTEX, particles.shape[0], range(particles.shape[0]))
    for p in range(Npart):
      out.InsertNextCell(vtk.VTK_LINE, Ntime, range(p, p + Npart*Ntime, Npart) )


  @staticmethod
  def request_information(source):
   pass


class ProgrammableDensityLoad(object):

  @staticmethod
  def main_script(source, density, extents=None, aname="default"):
    import vtk
    from vtk.util import numpy_support
  
    if len(density.shape) > 3:
      _, Nx, Ny, Nz = density.shape
    else:
      Nx, Ny, Nz = density.shape

    ido = source.GetOutput()
    ido.SetDimensions(Nx, Ny, Nz)
    if not extents is None:
      origin = extents[:6:2]
      spacing = (extents[1]-extents[0])/Nx, (extents[3]-extents[2])/Ny, (extents[5]-extents[4])/Nz
    else:
      origin = (-1, -1, -1)
      spacing = 2.0 / Nx, 2.0/Ny, 2.0/Nz
      
    ido.SetOrigin(*origin)
    ido.SetSpacing(*spacing)
    ido.SetExtent([0,Nx-1,0,Ny-1,0,Nz-1])
    if len(density.shape) > 3 and density.shape[0] == 3:
       N = Nx*Ny*Nz
       density = density.transpose().astype(np.float64).reshape((N,3))
       arr = numpy_support.numpy_to_vtk(density, deep=1)
    else:
       arr = numpy_support.numpy_to_vtk(density.transpose().astype(np.float64).ravel(), deep=1)
    arr.SetName(aname)
    ido.GetPointData().AddArray(arr)
    
  @staticmethod
  def request_information(source, density=None, dims=None):
    import vtk
    
    Nx = Ny = Nz = None
    if not density is None:
      Nx, Ny, Nz = density.shape
    elif not dims is None:
      Nx, Ny, Nz = dims
    else:
      raise ValueError("Need at least a density or dims")

    source.GetExecutive().GetOutputInformation(0).Set(vtk.vtkStreamingDemandDrivenPipeline.WHOLE_EXTENT(), 0, Nx-1, 0, Ny-1, 0, Nz-1)

  @staticmethod
  def prepare_timesteps_info(algorithm, timesteps):

    def SetOutputTimesteps(algorithm, timesteps):
      executive = algorithm.GetExecutive()
      outInfo = executive.GetOutputInformation(0)
      outInfo.Remove(executive.TIME_STEPS())
      for timestep in timesteps:
        outInfo.Append(executive.TIME_STEPS(), timestep)
      outInfo.Remove(executive.TIME_RANGE())
      outInfo.Append(executive.TIME_RANGE(), timesteps[0])
      outInfo.Append(executive.TIME_RANGE(), timesteps[-1])

    SetOutputTimesteps(algorithm, timesteps)

  @staticmethod
  @contextmanager
  def get_timestep(algorithm):

    def GetUpdateTimestep(algorithm):
      """Returns the requested time value, or None if not present"""
      executive = algorithm.GetExecutive()
      outInfo = executive.GetOutputInformation(0)
      if not outInfo.Has(executive.UPDATE_TIME_STEP()):
          return None
      return outInfo.Get(executive.UPDATE_TIME_STEP())

    # This is the requested time-step. This may not be exactly equal to the 
    # timesteps published in RequestInformation(). Your code must handle that
    # correctly
    req_time = GetUpdateTimestep(algorithm)

    output = algorithm.GetOutput()

    yield req_time
    
    # Now mark the timestep produced.
    output.GetInformation().Set(output.DATA_TIME_STEP(), req_time)

    
def load_borg(pdo, restart_name, mcmc_name, info=False, aname="BORG"):
  import h5py as h5
  
  with h5.File(restart_name) as f:
    N0 = f["/info/scalars/N0"][:]
    N1 = f["/info/scalars/N1"][:]
    N2 = f["/info/scalars/N2"][:]
    L0 = f["/info/scalars/L0"][:]
    L1 = f["/info/scalars/L1"][:]
    L2 = f["/info/scalars/L2"][:]
    c0 = f["/info/scalars/corner0"][:]
    c1 = f["/info/scalars/corner1"][:]
    c2 = f["/info/scalars/corner2"][:]

  if not info:
    with h5.File(mcmc_name) as f:
      d = f["/scalars/BORG_final_density"][:]+1
                
    ProgrammableDensityLoad.main_script(pdo, d, extents=[c0,c0+L0,c1,c1+L1,c2,c2+L2], aname=aname)
  else:
    ProgrammableDensityLoad.request_information(pdo, dims=[N0,N1,N2])
    
    
def load_borg_galaxies(pdo, restart_name, cid=0, info=False, aname="Galaxies"):
  import h5py as h5
  
  with h5.File(restart_name) as f:
    gals = f['/info/galaxy_catalog_%d/galaxies' % cid]
    ra = gals['phi'][:]
    dec = gals['theta'][:]
    r = gals['r'][:]

  if not info:
    x = r * np.cos(ra)*np.cos(dec)
    y = r * np.sin(ra)*np.cos(dec)
    z = r * np.sin(dec)
    parts = np.array([x,y,z]).transpose()
    ProgrammableParticleLoad.main_script(pdo, parts)
