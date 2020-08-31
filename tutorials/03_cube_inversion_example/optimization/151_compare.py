import numpy as np
import matplotlib.pyplot as plt
import os
from srfpython.HerrMet.files import ROOTNAME, HERRMETEXTRACTPDFMODELFILE
from srfpython.depthdisp.depthmodels import depthmodel_from_mod96


with np.load('../model/model.npz') as loader:
    x = loader['x']
    y = loader['y']
    z = loader['z']
    vstruth = loader['vs']


x0, nx, dx = x[0], len(x), x[1] - x[0]
y0, ny, dy = y[0], len(y), y[1] - y[0]
z0, nz, dz = z[0], len(z), z[1] - z[0]
# ========================== Display

vsinv = np.zeros_like(vstruth)
for iy in range(ny):
    for ix in range(nx):
        medfile = HERRMETEXTRACTPDFMODELFILE.format(
            rootname=ROOTNAME.format(node="node_iy{:03d}_ix{:03d}".format(iy, ix)),
            extract_mode="best",
            extract_limit=1000,
            extract_llkmin=0,
            extract_step=1,
            percentile=0.5)
        medfile = "../inversion/" +medfile
        assert os.path.isfile(medfile)

        dm = depthmodel_from_mod96(medfile)
        dm = dm.interp(ztop=z)
        vsinv[:, iy, ix] = dm.vs.values


xedges = x0 + np.arange(nx + 1) * dx - dx / 2
yedges = x0 + np.arange(ny + 1) * dy - dy / 2
zedges = z0 + np.arange(nz + 1) * dz - dz / 2; zedges[0] = 0.

for vs in vstruth, vsinv:
    # plt.figure()
    # for i in range(ny):
    #     for j in range(nx):
    #         plt.plot(z.flat[:], vs[:, i, j], 'r')
    # plt.gca().set_xlabel('z (km)')

    plt.figure()
    plt.colorbar(
        plt.pcolormesh(xedges, yedges, vs[0, ...],
                       vmin = 0.5, vmax = 3.5, cmap=plt.get_cmap('jet_r'))
        )
    plt.gca().set_xlabel('x (km)')
    plt.gca().set_ylabel('y (km)')

    plt.figure()
    plt.colorbar(
        plt.pcolormesh(xedges, zedges, vs[:, ny // 2, :],
                       vmin = 0.5, vmax = 3.5, cmap=plt.get_cmap('jet_r'))
        )
    plt.gca().set_xlabel('x (km)')
    plt.gca().set_ylabel('z (km)')
    plt.gca().invert_yaxis()



with np.load('mprior.npz') as loader:

    nz = loader["nz"]
    ny = loader["ny"]
    nx = loader["nx"]

    z = loader['ztop']
    dz = loader["ztop"][1] - loader["ztop"][0]
    dy = loader["y"][1] - loader["y"][0]
    dx = loader["x"][1] - loader["x"][0]

    vsopt = np.load('Msol.npy').reshape((nz, ny, nx))


xedges = x0 + np.arange(nx + 1) * dx - dx / 2
yedges = x0 + np.arange(ny + 1) * dy - dy / 2
zedges = z0 + np.arange(nz + 1) * dz - dz / 2; zedges[0] = 0.

vs = vsopt


plt.figure()
plt.colorbar(
    plt.pcolormesh(xedges, yedges, vs[0, ...],
                   vmin = 0.5, vmax = 3.5, cmap=plt.get_cmap('jet_r'))
    )
plt.gca().set_xlabel('x (km)')
plt.gca().set_ylabel('y (km)')

plt.figure()
plt.colorbar(
    plt.pcolormesh(xedges, zedges, vs[:, ny // 2, :],
                   vmin = 0.5, vmax = 3.5, cmap=plt.get_cmap('jet_r'))
)
plt.gca().set_xlabel('x (km)')
plt.gca().set_ylabel('z (km)')
plt.gca().invert_yaxis()

plt.ion()
plt.show()
input('pause')





