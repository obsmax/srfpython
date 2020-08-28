import numpy as np
from srfpython.depthdisp.depthmodels import depthmodel_from_arrays
import os
os.system('mkdir --parents ./model/nodes')

with np.load('./model/model.npz') as loader:
    x = loader['x']
    y = loader['y']
    z = loader['z']
    vs = loader['vs']
    vp = loader['vp']
    rh = loader['rh']

nx, ny, nz = len(x), len(y), len(z)

for iy in range(ny):
    for ix in range(nx):
        nodename = "node_iy%03d_ix%03d" % (iy, ix)
        filename = "./model/nodes/{}.mod96".format(nodename)

        dm = depthmodel_from_arrays(
            z=z,
            vp=vp[:, iy, ix],
            vs=vs[:, iy, ix],
            rh=rh[:, iy, ix])
        print(filename)
        dm.write96(filename=filename)
