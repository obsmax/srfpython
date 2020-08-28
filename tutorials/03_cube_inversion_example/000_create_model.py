from __future__ import print_function
from builtins import input
import matplotlib.pyplot as plt
import numpy as np
from srfpython.depthdisp.depthmodels import brocher2005
import os
os.system('mkdir --parents ./model')
"""
build a 3D vs model (vp, rh from Brocher 2005)
"""

x0 = 0.
y0 = 0.
z0 = 0.
nx = 15
ny = 10
nz = 10
dx = 2.4  # km
dy = 2.2  # km
dz = 0.33  # km

vs = np.zeros((nz, ny, nx), float)

z = (z0 + np.arange(nz) * dz)[:, np.newaxis, np.newaxis]
y = (y0 + np.arange(ny) * dy)[np.newaxis, :, np.newaxis]
x = (x0 + np.arange(nx) * dx)[np.newaxis, np.newaxis, :]

vs = 0.8 + 1.0 * z + 0.4 * \
     np.cos(2. * np.pi * x / 20.) * \
     np.cos(2. * np.pi * y / 20.) * \
     np.cos(2. * np.pi * z / 4.)

vp, rh = brocher2005(vs)
print('vsmin, vsmax', vs.min(), vs.max())


# ========================== Save
np.savez(
    './model/model.npz',
    x=x.flat[:],
    y=y.flat[:],
    z=z.flat[:],
    vs=vs,
    vp=vp,
    rh=rh)

# ========================== Display
xedges = x0 + np.arange(nx + 1) * dx - dx / 2
yedges = x0 + np.arange(ny + 1) * dy - dy / 2
zedges = z0 + np.arange(nz + 1) * dz - dz / 2; zedges[0] = 0.

plt.figure()
for i in range(ny):
    for j in range(nx):
        plt.plot(z.flat[:], vs[:, i, j], 'r')
        plt.plot(z.flat[:], vp[:, i, j], 'k')
        plt.plot(z.flat[:], rh[:, i, j], 'b')
plt.gca().set_xlabel('z (km)')


plt.figure()
plt.colorbar(
    plt.pcolormesh(xedges, yedges, vs[0, ...], cmap=plt.get_cmap('jet_r'))
    )
plt.gca().set_xlabel('x (km)')
plt.gca().set_ylabel('y (km)')

plt.figure()
plt.colorbar(
    plt.pcolormesh(xedges, zedges, vs[:, ny // 2, :], cmap=plt.get_cmap('jet_r'))
    )
plt.gca().set_xlabel('x (km)')
plt.gca().set_ylabel('z (km)')
plt.gca().invert_yaxis()

plt.ion()
plt.show()
input('pause')





