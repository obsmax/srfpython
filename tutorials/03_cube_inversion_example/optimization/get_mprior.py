import numpy as np
import matplotlib.pyplot as plt
#from sigy.utils.containers import Container
from srfpython.depthdisp.depthmodels import depthmodel1D, depthmodel_from_mod96
from srfpython.standalone.multipro8 import Job, MapSync
from srfpython.standalone.stdout import waitbarpipe

# ==================== read the parameter file
with open('optimize.param', 'r') as fid:
    for l in fid:
        l = l.split('\n')[0].split('#')[0].strip()
        if not len(l):
            continue
        if l.startswith('nx'):           nx = int(l.split()[1])
        if l.startswith('ny'):           ny = int(l.split()[1])
        if l.startswith('newnz'):     newnz = int(l.split()[1])
        if l.startswith('ndecim'):   ndecim = int(l.split()[1])

        if l.startswith('dx'):           dx = float(l.split()[1])
        if l.startswith('dy'):           dy = float(l.split()[1])
        if l.startswith('newdz'):     newdz = float(l.split()[1])
        if l.startswith('Lh'):           Lh = float(l.split()[1])
        if l.startswith('Lv'):           Lv = float(l.split()[1])
        if l.startswith('vsunc'):     vsunc = float(l.split()[1])
        if l.startswith('damping'): damping = float(l.split()[1])

        if l.startswith('extractfiles'): extractfiles = l.split()[1]
        if l.startswith('datafiles'):       datafiles = l.split()[1]


x = (np.arange(nx) * dx)[::ndecim]
y = (np.arange(ny) * dy)[::ndecim]
ixs = np.arange(nx)[::ndecim]  # indexs used to name the files
iys = np.arange(ny)[::ndecim]
ztop = np.arange(newnz) * newdz
# WARNING : I update nx, ny
nz, ny, nx = len(ztop), len(y), len(x)

# =============== load the results of the point wise inversion
# config = z, y, x
zmid = np.hstack((ztop[:-1] + 0.5 * (ztop[1:] - ztop[:-1]), ztop[-1] + ztop[1] - ztop[0]))
vs_cube = np.zeros((len(ztop), len(y), len(x)), float)
vsunc_cube = np.zeros_like(vs_cube)

for i, iy in enumerate(iys):
    for j, ix in enumerate(ixs):
        try:
            m96file = extractfiles.format(iy=iy, ix=ix)
            print('loading {}'.format(m96file))
            dm = depthmodel_from_mod96(m96file)
        except Exception as e:
            raise e

        vs_cube[:, i, j] = dm.vs.interp(z=zmid)
        vsunc_cube[:, i, j] = vsunc

# ===============
iz, jy, kx = np.mgrid[:nz, :ny, :nx]

ixsflat = ixs[kx].flat[:]
iysflat = iys[jy].flat[:]
xflat = x[kx].flat[:]
yflat = y[jy].flat[:]
ztopflat = ztop[iz].flat[:]
zmidflat = zmid[iz].flat[:]
Mprior = vs_cube   #.flat[:]
Munc = vsunc_cube  #.flat[:]

x_edges = np.hstack((x - 0.5 * (x[1] - x[0]), x[-1] + 0.5 * (x[1] - x[0])))
y_edges = np.hstack((y - 0.5 * (y[1] - y[0]), y[-1] + 0.5 * (y[1] - y[0])))
z_edges = np.hstack((ztop, ztop[-1] + 0.5 * (ztop[1] - ztop[0])))


# # ================== compute the parameterizers
# Mprior = Mprior.reshape((nz, ny, nx))
#
#
def get_parameterizers(verbose=True, **mapkwargs):
    """
    construct parameterizers needed to define the theory function in each node
    :param self:
    :return:
    """
    # write the header of the parameterizer for each node
    parameter_string_header = """
    #met NLAYER = {}
    #met TYPE = 'mZVSVPvsRHvp'
    #met VPvs = 'lambda VS: 0.9409+2.0947*VS-0.8206*VS**2+0.2683*VS**3-0.0251*VS**4'
    #met RHvp = 'lambda VP: 1.6612*VP-0.4721*VP**2+0.0671*VP**3-0.0043*VP**4+0.000106*VP**5'
    #fld KEY     VINF          VSUP
    #unt []      []            []
    #fmt %s      %f            %f
    """.format(len(ztop)).replace('    #', '#')

    for i in range(1, len(ztop)):
        # force VINF=VSUP => means lock the depth of the interfaces in the theory operator
        parameter_string_header += "-Z{} {} {}\n".format(i, -ztop[i], -ztop[i])  # add locked depth interfaces

    def job_generator():
        for iy in range(ny):
            for jx in range(nx):# order matters!!!!
                vs = Mprior[:, iy, jx]
                yield Job(iy, jx, ztop, vs)

    def job_handler(iy, jx, ztop, vs):
        parameterizer_string = parameter_string_header

        for i in range(len(ztop)):
            # SET VINF < VS extracted from pointwise inv < VSUP
            # such as parameterizer.MMEAN corresponds to the extracted vs
            parameterizer_string += "VS{} {} {}\n".format(i, vs[i] - 0.01, vs[i] + 0.01)
        # parameterizer = load_paramfile(parameter_string, verbose=False)[0]
        return iy, jx, parameterizer_string

    wb = None
    if verbose:
        wb = waitbarpipe('parameterizers')

    parameterizer_strings = []
    with MapSync(job_handler, job_generator(), **mapkwargs) as ma: # order matters!!!!
        for jobid, (iy, jx, parameter_string), _, _ in ma:
            parameterizer_strings.append(parameter_string)

            if verbose:
                wb.refresh(jobid / float(nx * ny))

    if verbose:
        wb.close()

    return np.asarray(parameterizer_strings, str)


parameterizer_strings = get_parameterizers(verbose=True)

np.savez(
    'mprior.npz',
    nz=nz, ny=ny, nx=nx,
    x=x, x_edges=x_edges, xflat=xflat, ixsflat=ixsflat,
    y=y, y_edges=y_edges, yflat=yflat, iysflat=iysflat,
    ztop=ztop, zmid=zmid, ztopflat=ztopflat, zmidflat=zmidflat,
    Mprior=Mprior.flat[:], Munc=Munc.flat[:],
    damping=damping, Lv=Lv, Lh=Lh,
    parameterizer_strings=parameterizer_strings)

# exit()
#
#
# plt.figure()
# iz = np.argmin(np.abs(zmid - 1.0))
# plt.subplot(221, aspect=1.0)
# plt.colorbar(
#     plt.pcolormesh(x, y, Mprior.reshape((nz, ny, nx))[iz, ...],
#                vmin=Mprior.min(), vmax=Mprior.max(),
#                cmap=plt.get_cmap('jet_r')))
#
# plt.subplot(222, aspect=1.0)
# plt.colorbar(
#     plt.pcolormesh(x, y, Munc.reshape((nz, ny, nx))[iz, ...],
#                vmin=Munc.min(), vmax=Munc.max(),
#                cmap=plt.get_cmap('nipy_spectral')))
#
# iy = np.argmin(np.abs(y - 10.0))
# plt.subplot(223)
# plt.colorbar(
#     plt.pcolormesh(x, zmid, Mprior.reshape((nz, ny, nx))[:, iy, :],
#                vmin=Mprior.min(), vmax=Mprior.max(),
#                cmap=plt.get_cmap('jet_r')))
# plt.gca().invert_yaxis()
#
# plt.subplot(224)
# plt.colorbar(
#     plt.pcolormesh(x, zmid, Munc.reshape((nz, ny, nx))[:, iy, :],
#                vmin=Munc.min(), vmax=Munc.max(),
#                cmap=plt.get_cmap('nipy_spectral')))
# plt.gca().invert_yaxis()
#
# plt.show()
