import struct
import numpy as np

def readGrafic(filename):
  """This function reads a grafic file.
  
  Arguments:
    filename (str): the path to the grafic file
    
  Returns:
    a tuple containing:
      * the array held in the grafic file
      * the size of the box
      * the scale factor
      * the mean matter density :math:`\Omega_\mathrm{m}`
      * the dark energy density :math:`\Omega_\Lambda`
      * the hubble constant, relative to 100 km/s/Mpc
      * xoffset
      * yoffset
      * zoffset
  """ 
  with open(filename, mode="rb") as f:
        p = struct.unpack("IIIIffffffffI", f.read(4*11 + 2*4))
        checkPoint0, Nx, Ny, Nz, delta, xoff, yoff, zoff, scalefac, omega0, omegalambda0, h, checkPoint1 = p
        if checkPoint0 != checkPoint1 or checkPoint0 != 4*11:
          raise ValueError("Invalid unformatted access")
          
        a = np.empty((Nx,Ny,Nz), dtype=np.float32)
        
        BoxSize = delta * Nx * h
        xoff *= h
        yoff *= h
        zoff *= h
        
        checkPoint = 4*Ny*Nz
        for i in range(Nx):
            checkPoint = struct.unpack("I", f.read(4))[0]
            if checkPoint != 4*Ny*Nz:
              raise ValueError("Invalid unformatted access")
              
            a[i, :, :] = np.fromfile(f, dtype=np.float32, count=Ny*Nz).reshape((Ny, Nz))

            checkPoint = struct.unpack("I", f.read(4))[0]
            if checkPoint != 4*Ny*Nz:
              raise ValueError("Invalid unformatted access")

  return a, BoxSize, scalefac, omega0, omegalambda0, h, xoff, yoff,zoff

def writeGrafic(filename, field, BoxSize, scalefac, **cosmo):

    with open(filename, mode="wb") as f:
        checkPoint = 4*11
        
        Nx,Ny,Nz = field.shape
        delta = BoxSize/Nx/cosmo['h']
        bad = 0.0
        
        f.write(struct.pack("IIIIffffffffI", checkPoint, 
                                          Nx, Ny, Nz, 
                                          delta, 
                                          bad, bad, bad, 
                                          scalefac,
                                          cosmo['omega_M_0'], cosmo['omega_lambda_0'], 100*cosmo['h'], checkPoint))
        checkPoint = 4*Ny*Nz
        field = field.reshape(field.shape, order='F')
        for i in range(Nx):
            f.write(struct.pack("I", checkPoint))
            f.write(field[i].astype(np.float32).tostring())
            f.write(struct.pack("I", checkPoint))
        
        
def writeWhitePhase(filename, field):

    with open(filename, mode="wb") as f:
        Nx,Ny,Nz = field.shape
        
        checkPoint = 4*4
        f.write(struct.pack("IIIIII", checkPoint, Nx, Ny, Nz, 0, checkPoint))
        
        field = field.reshape(field.shape, order='F')
        checkPoint = struct.pack("I", 4*Ny*Nz)
        for i in range(Nx):
            f.write(checkPoint)
            f.write(field[i].astype(np.float32).tostring())
            f.write(checkPoint)


def readWhitePhase(filename):
    with open(filename, mode="rb") as f:
        _, Nx, Ny, Nz, _, _ = struct.unpack("IIIIII", f.read(4*4+2*4))
        
        a = np.empty((Nx,Ny,Nz), dtype=np.float32)
        
        checkPoint_ref = 4*Ny*Nz
        
        for i in range(Nx):
            if struct.unpack("I", f.read(4))[0] != checkPoint_ref:
              raise ValueError("Invalid unformatted access")
              
            b = np.fromfile(f, dtype=np.float32, count=Ny*Nz).reshape((Ny, Nz))
            if i==0:
              print(b)
            a[i, : ,:] = b
            
            if struct.unpack("I", f.read(4))[0] != checkPoint_ref:
              raise ValueError("Invalid unformatted access")

    return a
