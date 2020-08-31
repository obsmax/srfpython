from __future__ import print_function
from builtins import input
import sys, glob, os
import warnings
import numpy as np
import scipy.sparse as sp
from scipy.sparse import linalg as splinalg
from srfpython.depthdisp.depthmodels import depthmodel_from_mod96
from srfpython.standalone.multipro8 import Job, MapSync
from srfpython.depthdisp.dispcurves import surf96reader
from srfpython.standalone.stdout import waitbarpipe
from srfpython.HerrMet.datacoders import Datacoder_log, makedatacoder
from srfpython.HerrMet.optimizetools import ForwardOperator, ModelSmoother, ModelCovarianceMatrix
from srfpython.utils import Timer
from srfpython.HerrMet.files import \
    HERRMETOPTIMIZEPARAMFILE, HERRMETOPTIMIZEPRIORFILE, \
    HERRMETOPTIMIZEDOBSFILE, HERRMETOPTIMIZEINVFILE

# ------------------------------ defaults
default_option = None

# ------------------------------ autorized_keys
authorized_keys = ["-h", "-help",
                   "-temp", "-mprior", "-dobs", "-inv", "-show"]

# ------------------------------ help messages
short_help = "--optimize   3D optimization plugin"

long_help = """\
--optimize           perform 3D linearized inversion starting from the 
                     output of a point-wise depth inversion
    -temp            create a template parameter file and exit (other options ignored)
    -mprior          collect the prior model from the point-wise depth inversion
                     load the files based on the content of the parameterfile 
                     created and customized after option -temp
    -dobs            collect the observed dispersion points to fit
    -inv             run the linearized inversion
    -show  s f       show a slice accross the prior (from point wise depth inv), 
                     starting and final models
                     provide the slice direction z y or x and the slice value in km
    -h, -help        display the help message for this plugin 
""".format(default_option=default_option)

# ------------------------------ example usage
example = """\
## OPTIMIZE
# ... 

HerrMet --optimize 
"""


def check_keys(argv):
    for k in argv.keys():
        if k in ['main', "_keyorder"]:
            continue  # private keys

        if k not in authorized_keys:
            print('option {} is not recognized'.format(k))  # only for the default plugin

            raise Exception('option {} is not recognized'.format(k))  # please use this line in other plugins


OPTIMIZEPARAMTXT = """
# Parameter file for the optimization plugin
# ========== recall the grid parameters
nx           {nx:d}
ny           {ny:d}
dx           {dx:f}    # km
dy           {dy:f}    # km

# ========== path to the point-wise inversion data and results
# formattable strings including reference to iy and ix
targetfiles  {targetfiles:s}    
extractfiles {extractfiles:s}   

# ========== optimization parameters
ndecim       {ndecim:d}  # downsamp the grid
newnz        {newnz:f}   # new vertical grid (needed to interpolate the extracted files)
newdz        {newdz:f}   # new vertical grid (needed to interpolate the extracted files)
Lh           {Lh:f}   # horizontal smoothing distance km
Lv           {Lv:f}   # vertical smoothing distance km
vsunc        {vsunc:f}   # km/s
damping      {damping:f}   # weight of the regularization
"""


def load_optimize_parameters():
    with open(HERRMETOPTIMIZEPARAMFILE, 'r') as fid:
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

    return nx, ny, newnz, dx, dy, newdz, ndecim, Lh, Lv, vsunc, damping, extractfiles, datafiles


def add_to_npz(npzfile, allow_pickle=False, **kwargs):
    """
    add onei or more field(s) to an existing npzfile (implies loading the existing fields)
    :param npzfile: name of the file to load and write
    :param kwargs: fields to add or update
    """
    if os.path.isfile(npzfile):

        # load if not in kwargs
        keys = kwargs.keys()
        with np.load(npzfile) as loader:
            loaded = {f: loader[f] for f in loader.files if f not in keys}

        # concatenate
        d = dict(loaded, **kwargs)

    else:
        d = kwargs

    # save
    from numpy.lib.npyio import _savez
    _savez(npzfile, args=(), kwds=d, compress=True, allow_pickle=allow_pickle)


def optimize(argv, verbose, mapkwargs):

    if '-h' in argv.keys() or "-help" in argv.keys():
        print(long_help)
        return

    check_keys(argv)

    # ===============================
    if "-temp" in argv.keys():
          template = OPTIMIZEPARAMTXT.format(
              nx=1, ny=1, newnz=1,
              dx=1., dy=1.0, newdz=1.0,
              targetfiles="./path/to/inversion/_HerrMet_node_iy{iy:03d}_ix{ix:03d}/_HerrMet.target",
              extractfiles="./path/to/inversion/_HerrMet_node_iy{iy:03d}_ix{ix:03d}/_HerrMet.best_1000_0_1.p0.50.mod96",
              ndecim=1, Lh=1.0, Lv=0.5, vsunc=0.1, damping=1.0)
          print(template)

          if os.path.isfile(HERRMETOPTIMIZEPARAMFILE):
              raise IOError(HERRMETOPTIMIZEPARAMFILE, ' exists already')

          with open(HERRMETOPTIMIZEPARAMFILE, 'w') as fid:
              fid.write(template)

          print('please customize {} and rerun this plugin'.format(HERRMETOPTIMIZEPARAMFILE))
          sys.exit(0)

    # ===============================
    if "-mprior" in argv.keys():

        if not os.path.isfile(HERRMETOPTIMIZEPARAMFILE):
            raise IOError(HERRMETOPTIMIZEPARAMFILE)

        # ==================== read the parameter file
        nx, ny, newnz, \
            dx, dy, newdz, ndecim, \
            Lh, Lv, vsunc, damping, \
            extractfiles, datafiles = \
            load_optimize_parameters()

        x = (np.arange(nx) * dx)[::ndecim]
        y = (np.arange(ny) * dy)[::ndecim]
        ixs = np.arange(nx)[::ndecim]  # indexs used to name the files
        iys = np.arange(ny)[::ndecim]
        ztop = np.arange(newnz) * newdz
        # WARNING : I update nx, ny
        del dx, dy  # to avoid confusion
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
        Mprior = vs_cube  # .flat[:]
        Munc = vsunc_cube  # .flat[:]

        xedges = np.hstack((x - 0.5 * (x[1] - x[0]), x[-1] + 0.5 * (x[1] - x[0])))
        yedges = np.hstack((y - 0.5 * (y[1] - y[0]), y[-1] + 0.5 * (y[1] - y[0])))
        zedges = np.hstack((ztop, ztop[-1] + 0.5 * (ztop[1] - ztop[0])))

        # ================== compute the parameterizers
        """
        construct parameterizers needed to define the theory function in each node
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
                for jx in range(nx):  # order matters!!!!
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
        with MapSync(job_handler, job_generator(), **mapkwargs) as ma:  # order matters!!!!
            for jobid, (iy, jx, parameter_string), _, _ in ma:
                parameterizer_strings.append(parameter_string)

                if verbose:
                    wb.refresh(jobid / float(nx * ny))

        if verbose:
            wb.close()

        parameterizer_strings = np.asarray(parameterizer_strings, str)

        np.savez(
            HERRMETOPTIMIZEPRIORFILE,
            nz=nz, ny=ny, nx=nx,
            x=x, xedges=xedges, xflat=xflat, ixs=ixs, ixsflat=ixsflat,
            y=y, yedges=yedges, yflat=yflat, iys=iys, iysflat=iysflat,
            ztop=ztop, zmid=zmid, zedges=zedges, ztopflat=ztopflat, zmidflat=zmidflat,
            Mprior=Mprior, Munc=Munc,
            damping=damping, Lv=Lv, Lh=Lh,
            parameterizer_strings=parameterizer_strings)

    # ===============================
    if "-dobs" in argv.keys():

        # ==================== read the parameter file
        nx, ny, newnz, \
            dx, dy, newdz, ndecim, \
            Lh, Lv, vsunc, damping, \
            extractfiles, datafiles = \
            load_optimize_parameters()

        # x = (np.arange(nx) * dx)[::ndecim]
        # y = (np.arange(ny) * dy)[::ndecim]
        ixs = np.arange(nx)[::ndecim]  # indexs used to name the files
        iys = np.arange(ny)[::ndecim]

        datacoder_strings = []
        for i, iy in enumerate(iys):
            for j, ix in enumerate(ixs):  # order matters!!!!

                try:
                    datafile = datafiles.format(iy=iy, ix=ix)
                    if verbose:
                        print('loading ', datafile)
                    s96 = surf96reader(filename=datafile)

                except Exception as e:
                    raise e

                datacoder_string = str(s96)
                datacoder_strings.append(datacoder_string)

        datacoder_strings = np.asarray(datacoder_strings, str)
        datacoders = [makedatacoder(datacoder_string, which=Datacoder_log) for datacoder_string in
                      datacoder_strings]

        Nobs = []
        Waves = []
        Types = []
        Modes = []
        Freqs = []
        Dobs = []
        Dunc = []
        for datacoder in datacoders:
            _dobs, _CDinv = datacoder.target()
            Nobs.append(len(datacoder.waves))
            Waves.append(datacoder.waves)
            Types.append(datacoder.types)
            Modes.append(datacoder.modes)
            Freqs.append(datacoder.freqs)
            Dobs.append(_dobs)
            Dunc.append(_CDinv ** -0.5)

        Nobs = np.array(Nobs)
        Waves = np.hstack(Waves)
        Types = np.hstack(Types)
        Modes = np.hstack(Modes)
        Freqs = np.hstack(Freqs)
        Dobs = np.hstack(Dobs)
        Dunc = np.hstack(Dunc)

        np.savez(HERRMETOPTIMIZEDOBSFILE,
                 datacoder_strings=datacoder_strings,
                 Nobs=Nobs, Waves=Waves,
                 Types=Types, Modes=Modes,
                 Freqs=Freqs, Dobs=Dobs, Dunc=Dunc)

    # ===============================
    if "-inv" in argv.keys():

        # ========== Load inputs
        with np.load(HERRMETOPTIMIZEPRIORFILE) as loader:
            x = loader['x']
            y = loader['y']
            # xedges = loader['xedges']
            # yedges = loader['yedges']
            # zedges = loader['zedges']
            # ixs = loader['ixs']
            # iys = loader['iys']
            zmid = loader['zmid']
            xflat = loader['xflat']
            yflat = loader['yflat']
            zflat = zmidflat = loader['zmidflat']
            Lh = loader['Lh']
            Lv = loader['Lv']
            Mprior = loader['Mprior'].flat[:]
            Munc = loader['Munc'].flat[:]
            damping = loader['damping']
            parameterizer_strings = loader['parameterizer_strings']

        with np.load(HERRMETOPTIMIZEDOBSFILE) as loader:
            Nobs = loader['Nobs']
            # Waves = loader['Waves']
            # Types = loader['Types']
            # Modes = loader['Modes']
            # Freqs = loader['Freqs']
            Dobs = loader['Dobs']
            Dunc = loader['Dunc']
            datacoder_strings = loader['datacoder_strings']

        nx, ny, nz = len(x), len(y), len(zmid)
        nobs = len(Dobs)

        # save the grid
        # add_to_npz(HERRMETOPTIMIZEINVFILE,
        #            # 3D grid needed for display
        #            nx=nx, ny=ny, nz=nz,
        #            xedges=xedges, yedges=yedges, zedges=zedges,  # edge grid points (for pcolormesh)
        #            x=x, y=y, zmid=zmid,  # mid grid points
        #            ixs=ixs, iys=iys,  # indexs of the point wise inversion used
        #            Nobs=Nobs)   # number of disp point per node (array) flatten in the horizontal plane

        Munc *= np.sqrt(float(len(Mprior)) / float(len(Dobs))) / damping  # equilibrate the costs, apply damping here

        CD = sp.diags(Dunc ** 2., format="csc", shape=(len(Dobs), len(Dobs)))
        CDinvdiag = Dunc ** -2.

        smoother = ModelSmoother(
            x=x, y=y, zmid=zmid,
            xflat=xflat, yflat=yflat, zmidflat=zflat,
            Lh=Lh, Lv=Lv)

        CM = ModelCovarianceMatrix(
            Munc=Munc,
            x=x, y=y, zmid=zmid,
            xflat=xflat, yflat=yflat, zmidflat=zflat,
            Lh=Lh, Lv=Lv)

        g = ForwardOperator(
            parameterizer_strings=parameterizer_strings,
            datacoder_strings=datacoder_strings,
            nx=nx, ny=ny, nz=nz, Nobs=Nobs,
            **mapkwargs)

        Dprior = g(Mprior)
        add_to_npz(HERRMETOPTIMIZEINVFILE,
                   Dprior=Dprior)

        M0 = smoother.dot(Mprior, trunc=4.)

        D0 = g(M0)
        add_to_npz(HERRMETOPTIMIZEINVFILE,
                   M0=M0.reshape((nz, ny, nx)),
                   D0=D0)

        MI = M0.copy()
        DI = D0.copy()

        with Timer('CM_sparse'):
            CM_sparse = CM.sparse(trunc=4.)
        with Timer('CM_lu'):
            CM_sparse_lu = splinalg.splu(CM_sparse)

        G = None
        data_costs = []
        model_costs = []
        for i in range(20):
            if G is None or True:
                G = g.frechet_derivatives(MI)
                with Timer('G * CM * GT + CD'):
                    Ksparse = G * CM_sparse * G.T + CD

                with Timer('splu(G * CM * GT + CD)'):
                    Klu = splinalg.splu(Ksparse)

            XI = Dobs - DI + G * (MI - M0)
            DM = CM.dot(G.T * Klu.solve(XI))
            DM = (M0 - MI + DM)

            if np.abs(DM).max() > 0.25:
                warnings.warn(str(np.abs(DM).max()))
                DM *= 0.25 / np.abs(DM).max()

            print("max correction:", np.abs(DM).max())
            MI = MI + DM
            DI = g(MI)

            # ============== costs
            chi2_data = 0.5 * (CDinvdiag * (Dobs - DI) ** 2.0).sum()
            chi2_model = 0.5 * (CM_sparse_lu.solve(MI - M0) * (MI - M0)).sum()  # fails if I use Mprior instead of M0
            if chi2_model < 0 or chi2_data < 0:
                raise ValueError(chi2_model, chi2_data)  # may happen after numerical instability

            print('data cost: ', chi2_data)
            print('model cost: ', chi2_model)
            print('total cost: ', chi2_data + chi2_model)

            data_costs.append(chi2_data)
            model_costs.append(chi2_model)

            add_to_npz(HERRMETOPTIMIZEINVFILE,
                data_costs=np.asarray(data_costs),
                model_costs=np.asarray(model_costs),
                total_costs=np.asarray(model_costs) + np.asarray(data_costs),
                Msol=MI.reshape((nz, ny, nx)),
                Dsol=DI)

            if np.abs(DM).max() <= 0.01:
                print('convergence achieved')
                break

    # ===============================
    if "-show" in argv.keys():
        import matplotlib.pyplot as plt

        with np.load(HERRMETOPTIMIZEPRIORFILE) as loader:
            x = loader['x']
            y = loader['y']
            zmid = loader['zmid']
            xedges = loader['xedges']
            yedges = loader['yedges']
            zedges = loader['zedges']
            vs_prior = loader['Mprior']

        with np.load(HERRMETOPTIMIZEINVFILE) as loader:
            vs_0 = loader['M0']
            vs_sol = loader['Msol']

        if argv['-show'][0] == "z":
            zslice = argv['-show'][1]
            iz = np.argmin(abs(zmid - zslice))

            fig = plt.figure(figsize=(12, 4))
            ax1 = plt.subplot(131, title="$ Vs^{prior} $",
                              aspect=1.0, xlabel="x(km)", ylabel="y(km)")
            ax2 = plt.subplot(132, title="$ Vs^{0} $",
                              aspect=1.0, xlabel="x(km)", ylabel="y(km)",
                              sharex=ax1, sharey=ax1)
            ax3 = plt.subplot(133, title="$ Vs^{sol} $",
                              aspect=1.0, xlabel="x(km)", ylabel="y(km)",
                              sharex=ax1, sharey=ax1)
            fig.suptitle('depth : {:.2f}km'.format(zmid[iz]))
            cax = fig.add_axes((0.95, 0.3, 0.012, 0.3))

            vsmin = min([vs[iz, ...].min() for vs in [vs_prior, vs_0, vs_sol]])
            vsmax = min([vs[iz, ...].max() for vs in [vs_prior, vs_0, vs_sol]])
            for ax, vs in zip([ax1, ax2, ax3], [vs_prior, vs_0, vs_sol]):
                coll = ax.pcolormesh(xedges, yedges, vs[iz, ...],
                              vmin=vsmin, vmax=vsmax,
                              cmap=plt.get_cmap('jet_r'))

            plt.colorbar(coll, cax=cax)
            plt.ion()
            plt.show()
            input('pause')

        elif argv['-show'][0] == "x":
            xslice = argv['-show'][1]
            ix = np.argmin(abs(x - xslice))

            fig = plt.figure(figsize=(12, 4))
            ax1 = plt.subplot(131, title="$ Vs^{prior} $",
                              xlabel="y(km)", ylabel="z(km)")
            ax2 = plt.subplot(132, title="$ Vs^{0} $",
                              xlabel="y(km)", ylabel="z(km)",
                              sharex=ax1, sharey=ax1)
            ax3 = plt.subplot(133, title="$ Vs^{sol} $",
                              xlabel="y(km)", ylabel="z(km)",
                              sharex=ax1, sharey=ax1)
            for ax in ax1, ax2, ax3:
                ax.invert_yaxis()
            fig.suptitle('x = {:.2f}km'.format(x[ix]))
            cax = fig.add_axes((0.95, 0.3, 0.012, 0.3))

            vsmin = min([vs[..., ix].min() for vs in [vs_prior, vs_0, vs_sol]])
            vsmax = min([vs[..., ix].max() for vs in [vs_prior, vs_0, vs_sol]])
            for ax, vs in zip([ax1, ax2, ax3], [vs_prior, vs_0, vs_sol]):
                coll = ax.pcolormesh(yedges, zedges, vs[..., ix],
                                     vmin=vsmin, vmax=vsmax,
                                     cmap=plt.get_cmap('jet_r'))

            plt.colorbar(coll, cax=cax)
            plt.ion()
            plt.show()
            input('pause')

        elif argv['-show'][0] == "y":
            yslice = argv['-show'][1]
            iy = np.argmin(abs(y - yslice))

            fig = plt.figure(figsize=(12, 4))
            ax1 = plt.subplot(131, title="$ Vs^{prior} $",
                              xlabel="y(km)", ylabel="z(km)")
            ax2 = plt.subplot(132, title="$ Vs^{0} $",
                              xlabel="y(km)", ylabel="z(km)",
                              sharex=ax1, sharey=ax1)
            ax3 = plt.subplot(133, title="$ Vs^{sol} $",
                              xlabel="y(km)", ylabel="z(km)",
                              sharex=ax1, sharey=ax1)
            for ax in ax1, ax2, ax3:
                ax.invert_yaxis()
            fig.suptitle('y = {:.2f}km'.format(y[iy]))
            cax = fig.add_axes((0.95, 0.3, 0.012, 0.3))

            vsmin = min([vs[:, iy, :].min() for vs in [vs_prior, vs_0, vs_sol]])
            vsmax = min([vs[:, iy, :].max() for vs in [vs_prior, vs_0, vs_sol]])
            for ax, vs in zip([ax1, ax2, ax3], [vs_prior, vs_0, vs_sol]):
                coll = ax.pcolormesh(xedges, zedges, vs[:, iy, :],
                                     vmin=vsmin, vmax=vsmax,
                                     cmap=plt.get_cmap('jet_r'))

            plt.colorbar(coll, cax=cax)
            plt.ion()
            plt.show()
            input('pause')
