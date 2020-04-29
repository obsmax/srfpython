from __future__ import print_function
from srfpython import *
import os

os.system('mkdir --parents models')
os.system('trash models/model[0-9][0-9][0-9].mod96')

# -----------------------
# create a 1-D depth model
# -----------------------
# top layers array first at 0, positive, growing, km
ztop = np.linspace(0., 2.8, 50)

# vs in km/s
vs = (3.5 - .86) / (ztop[-1] - ztop[0]) * (ztop - 0.) + .86 + \
     -.5 * np.exp(-ztop / .4) + \
     .08 * np.cos(2. * np.pi * ztop / .5) + \
     .2 * np.sin(2. * np.pi * ztop / 1.) + \
     .1 * np.cos(2. * np.pi * ztop / 2.) + \
     .15 * np.cos(2. * np.pi * ztop / 3.)

vp, rh = brocher2005(vs)
dm = depthmodel_from_arrays(ztop, vp, vs, rh)

# -----------------------
# nodes in the 2D plan
# -----------------------
longitudes = 5. + np.linspace(-0.5, 0.5, 6)
latitudes = 40. + np.linspace(-0.5, 0.5, 5)
longitudes, latitudes = np.meshgrid(longitudes, latitudes)
longitudes += np.random.randn(*longitudes.shape) * 0.01
latitudes += np.random.randn(*latitudes.shape) * 0.01

longitudes = longitudes.flat[:]
latitudes = latitudes.flat[:]
depthmodels = [dm.copy() for _ in range(len(longitudes))]

# -----------------------
# create 3-D depth model by perturbing the depth profile at each point
# -----------------------
lon_anom     = [ 4.75,   5.,    5.25]
lat_anom     = [ 39.75,  40.,   40.25]
depth_anom   = [ 0.5,    1.5,   2.5]
rel_amp_anom = [+0.40,  -0.30, -0.20]
sigma_lons   = [ 0.4,    0.4,   0.4]
sigma_lats   = [ 0.4,    0.4,   0.4]
sigma_depth  = [ 0.5,    0.5,   0.5]

for lon, lat, dm in zip(longitudes, latitudes, depthmodels):
    for lon0, lat0, z0, amp, slo, sla, sz in zip(lon_anom, lat_anom, depth_anom, rel_amp_anom, sigma_lons, sigma_lats, sigma_depth):
        dm.vs.values *= (1.0 + amp * np.exp(-.5 * ((lon - lon0) / slo) ** 2.0)
                                   * np.exp(-.5 * ((lat - lat0) / sla) ** 2.0)
                                   * np.exp(-.5 * ((ztop - z0) / sz) ** 2.0))
    dm.vp.values, dm.rh.values = brocher2005(dm.vs.values)

fig1 = plt.figure(figsize=(4, 5))
for dm in depthmodels:
    dm.show(plt.gca())


fig2 = plt.figure(figsize=(4 * len(depth_anom), 3))
for nz in range(len(depth_anom)):
    vs_values = np.array([dm.vs.interp(depth_anom[nz]) for dm in depthmodels], float)
    plt.subplot(1, len(depth_anom), nz + 1, title=str(depth_anom[nz]) + "km")
    plt.tripcolor(
        longitudes, latitudes, vs_values,
        vmin=vs_values.min(), vmax=vs_values.max(),
        cmap=plt.get_cmap("jet_r"))

    cb = plt.scatter(
        longitudes,
        latitudes,
        50.,
        vs_values,
        vmin=vs_values.min(), vmax=vs_values.max(),
        cmap=plt.get_cmap("jet_r"))
    plt.colorbar(cb)

fig3 = plt.figure(figsize=(8, 4))
fig3.subplots_adjust(bottom=0.2)
z_interp = np.linspace(0., 3.5, 100)
ax1 = fig3.add_subplot(111)
E, S = [],  []
for n, dm in enumerate(depthmodels):
    vs = dm.vs.interp(z_interp)
    E.append(vs)

cb = ax1.pcolormesh(
    np.arange(len(depthmodels) + 1) - 0.5,
    np.hstack((z_interp, z_interp[-1] + 0.1)),
    np.array(E).T,
    vmin=0.2, vmax=3.5, cmap=plt.get_cmap('nipy_spectral_r'))
plt.colorbar(cb)
ax1.invert_yaxis()
ax1.set_xticks(range(len(depthmodels)))
ax1.set_xticklabels(['node{:03d}'.format(n) for n in range(len(depthmodels))],
                    verticalalignment="top",
                    horizontalalignment='left',
                    rotation=-45.)


# -----------------------
# save the models
# -----------------------

for nmodel, dm in enumerate(depthmodels):
    f = "models/node{:03d}.mod96".format(nmodel)
    print(f)
    dm.write96(f)

with open('nodes.txt', 'w') as fid:
    fid.write('#longitude latitude node\n')
    for n, (lon, lat) in enumerate(zip(longitudes, latitudes)):
        fid.write('{} {} node{:03d}\n'.format(lon, lat, n))

fig1.savefig('depth_models.png')
fig2.savefig('vs_slices.png')
fig3.savefig('profile.png')