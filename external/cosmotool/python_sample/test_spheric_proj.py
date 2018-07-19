import cosmotool as ct
import numpy as np
import healpy as hp

d = np.zeros((64,64,64), ct.DTYPE)

d[32,32,32] = 1

ii=np.arange(256)*64/256.-32
xx = ii[:,None,None].repeat(256,axis=1).repeat(256,axis=2).reshape(256**3)
yy = ii[None,:,None].repeat(256,axis=0).repeat(256,axis=2).reshape(256**3)
zz = ii[None,None,:].repeat(256,axis=0).repeat(256,axis=1).reshape(256**3)

d_high = ct.interp3d(xx, yy, zz, d, 64, periodic=True)
d_high = d_high.reshape((256,256,256))

proj0 = ct.spherical_projection(64, d, 0, 20, integrator_id=0, shifter=np.array([0.5,0.5,0.5]))
proj1 = ct.spherical_projection(64, d, 0, 20, integrator_id=1)

proj0_high = ct.spherical_projection(256, d_high, 0, 30, integrator_id=0, shifter=np.array([0.5,0.5,0.5]))
