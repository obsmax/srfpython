from srfpython import *
import time
input = raw_input

A = np.loadtxt('nodes.txt', comments="#", dtype=str)
lons = np.asarray(A[:, 0], float)
lats = np.asarray(A[:, 1], float)
nodes = A[:, 2]

expected_models = []
solution_models = []
for lon, lat, node in zip(lons, lats, nodes):
    expected_model = "./models/{}.mod96".format(node)
    solution_model = "./inversion/_HerrMet_{}/_HerrMet.best_1000_0_1.p0.50.mod96".format(node)
    if not os.path.isfile(expected_model) or not os.path.isfile(solution_model):
        raise IOError(expected_model, solution_model)

    dm_expected = depthmodel_from_mod96(expected_model)
    dm_solution = depthmodel_from_mod96(solution_model)

    expected_models.append(dm_expected)
    solution_models.append(dm_solution)


fig = plt.figure(figsize=(6, 6))
fig.subplots_adjust(bottom=0.2)
z_interp = np.linspace(0., 3.5, 100)
ax1 = fig.add_subplot(211, title="expected vs (km/s)")
ax2 = fig.add_subplot(212, title="solution vs (km/s)", sharex=ax1, sharey=ax1)
E, S = [],  []
for n, (dm_expected, dm_solution) in enumerate(zip(expected_models, solution_models)):

    vsexpected = dm_expected.vs.interp(z_interp)
    vssolution = dm_solution.vs.interp(z_interp)
    E.append(vsexpected)
    S.append(vssolution)

vmin=0.2
vmax=3.5
cmap=plt.get_cmap('nipy_spectral_r')

ax1.pcolormesh(
    np.arange(len(expected_models) + 1) - 0.5,
    np.hstack((z_interp, z_interp[-1] + 0.1)),
    np.array(E).T,
    vmin=vmin, vmax=vmax, cmap=cmap)

cb = ax2.pcolormesh(
    np.arange(len(solution_models) + 1) - 0.5,
    np.hstack((z_interp, z_interp[-1] + 0.1)),
    np.array(S).T,
    vmin=vmin, vmax=vmax, cmap=cmap)

ax1.invert_yaxis()
plt.setp(ax1.get_xticklabels(), visible=False)
ax2.set_xticks(range(len(expected_models)))
ax2.set_xticklabels(['node{:03d}'.format(n) for n in range(len(expected_models))],
                    verticalalignment="top",
                    horizontalalignment='left',
                    rotation=-45.)
cax = fig.add_axes((0.92, 0.3, 0.012, 0.4))
plt.colorbar(cb, cax=cax)
fig.savefig('inverted_expected_profiles.png')


fig = plt.figure(figsize=(8, 8))
for n, z0 in enumerate([.5, 1.5, 2.5]):
    vs_expected = np.asarray([dm.vs.interp(z0) for dm in expected_models])
    vs_solution = np.asarray([dm.vs.interp(z0) for dm in solution_models])

    ax = fig.add_subplot(3, 2, 1 + 2 * n, title="expected vs at depth {}km".format(z0))
    bx = fig.add_subplot(3, 2, 2 + 2 * n, title="solution vs at depth {}km".format(z0))

    ax.tripcolor(
        lons, lats, vs_expected,
        vmin=vs_expected.min(),
        vmax=vs_expected.max(),
        cmap=plt.get_cmap('jet_r'))

    ax.scatter(lons, lats, 20, vs_expected,
               vmin=vs_expected.min(),
               vmax=vs_expected.max(),
               cmap=plt.get_cmap('jet_r'))

    bx.tripcolor(
        lons, lats, vs_solution,
        vmin=vs_expected.min(),
        vmax=vs_expected.max(),
        cmap=plt.get_cmap('jet_r'))

    coll = bx.scatter(lons, lats, 20, vs_solution,
               vmin=vs_expected.min(),
               vmax=vs_expected.max(),
               cmap=plt.get_cmap('jet_r'))

    plt.colorbar(coll, ax=bx)

fig.savefig('inverted_expected_slices.png')
plt.ion()
plt.show()
input('pause : press enter')

# plt.figure(figsize=(4, 6))
# for dm_expected, dm_solution in zip(expected_models, solution_models):
#     plt.gca().cla()
#     dm_expected.show(plt.gca(), linewidth=2)
#     dm_solution.show(plt.gca(), linewidth=1, linestyle="--")
#
#     plt.gca().set_title(node)
#     plt.ion()
#     plt.gcf().show()
#     plt.gcf().canvas.draw()
#     time.sleep(0.1)
#
# plt.show()