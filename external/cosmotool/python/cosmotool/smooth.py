from .config import install_prefix
import subprocess
import os
try:
  from tempfile import TemporaryDirectory
except:
  from backports.tempfile import TemporaryDirectory

import h5py as h5
import numpy as np
import weakref

def smooth_particle_density(
           position,
           velocities=None,
           radius=1e6,
           boxsize=None, 
           resolution=128,
           center=None, tmpprefix=None ):
  """Use adaptive smoothing to produce density and momentum fields.
  The algorithm is originally described in [1].

  Parameters:
  position   : numpy array NxQ
               the particle positions
               if Q==3, only positions. Q==6 means full space phase
  velocities : Optional numpy array Nx3. 
               It is only optional if the above Q is 6.
  radius     : float
               Maximum radius to which we need to compute fields
  boxsize    : float
               Size of the box for the generated fields
  resolution : int
               Resolution of the output boxes
  center     : list of 3 floats
               Center of the new box. It depends on the convention 
               for particles. If those are between [0, L], then [0,0,0]
               is correct. If those are [-L/2,L/2] then you should set
               [L/2,L/2,L/2].
  tmpprefix  : string
               prefix of the temporary directory that will be used.
               It needs to have a lot of space available. By default
               '/tmp/ will be typically used.


  Returns
  -------
   dictionnary
   The dict has two entries: 'rho' for the density, and 'p' for the momenta.
   Once the dictionary is garbage collected all temporary files and directories
   will be cleared automatically.


  Raises
  ------
   ValueError
     if arguments are invalid


  .. [1] S. Colombi, M. Chodorowski, 
         "Cosmic velocity-gravity in redshift space", MNRAS, 2007, 375, 1
  """
  if len(position.shape) != 2:
     raise ValueError("Invalid position array shape")

  if velocities is None:
    if position.shape[1] != 6:
      raise ValueError("Position must be phase space if no velocities are given")


  if boxsize is None:
    raise ValueError("Need a boxsize")

  cx,cy,cz=center
  tmpdir = TemporaryDirectory(prefix=tmpprefix) 
  h5_file = os.path.join(tmpdir.name, 'particles.h5')
  with h5.File(h5_file, mode="w") as f:
    data = f.create_dataset('particles', shape=(position.shape[0],7), dtype=np.float32)
    data[:,:3] = position[:,:3]
    if velocities is not None: 
      data[:,3:6] = velocities[:,:3]
    else:
      data[:,3:6] = position[:,3:]
    data[:,6] = 1


  ret =  \
     subprocess.run([
         os.path.join(install_prefix,'bin','simple3DFilter'),
         h5_file,
         str(radius),
         str(boxsize),
         str(resolution),
         str(cx), str(cy), str(cz)
     ], cwd=tmpdir.name)

  
  f0 = h5.File(os.path.join(tmpdir.name,'fields.h5'), mode="r")
  def cleanup_f0():
    f0.close()
    tmpdir.cleanup()

  class Dict(dict):
    pass

  t = Dict(rho=f0['density'], p=[f0['p0'], f0['p1'], f0['p2']])
  t._tmpdir_=tmpdir
  weakref.finalize(t, cleanup_f0)
  return t
