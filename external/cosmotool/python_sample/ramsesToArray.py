import numpy as np
import sys
import h5py as h5
import cosmotool as ct

s = ct.loadRamsesAll(sys.argv[1], int(sys.argv[2]), doublePrecision=True, loadVelocity=True)


q = [p for p in s.getPositions()]
q += [p for p in s.getVelocities()]
q += [np.ones(q[0].size,dtype=q[0].dtype)]
q = np.array(q)

with h5.File("particles.h5", mode="w") as f:
  f.create_dataset("particles", data=q.transpose())
