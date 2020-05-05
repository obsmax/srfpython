from __future__ import print_function
from builtins import input
import sys, glob, os
import warnings
import numpy as np
from scipy.sparse import spmatrix, diags, csr_matrix, csc_matrix, \
    save_npz as save_sparse_npz, load_npz as load_sparse_npz
from scipy.sparse.linalg import inv as sparse_inv
import matplotlib.pyplot as plt
# from srfpython.HerrMet.files import *
from srfpython.HerrMet.nodefile import NodeFile
from srfpython.depthdisp.depthmodels import depthmodel_from_mod96, depthmodel1D, depthmodel_from_arrays
from srfpython.HerrMet.paramfile import load_paramfile
from srfpython.HerrMet.datacoders import makedatacoder, Datacoder_log
from srfpython.HerrMet.theory import Theory
from srfpython.standalone.multipro8 import Job, MapSync


WORKINGDIR = "."
NODEFILELOCAL = os.path.join("{workingdir}", "_HerrMet_.nodes")
MPRIORFILE = os.path.join("{workingdir}", "_HerrMet_Mprior.npy")

MFILE = os.path.join("{workingdir}", "_HerrMet_M{niter:03d}.npy")
MFILES = os.path.join("{workingdir}", "_HerrMet_M[0-9][0-9][0-9].npy")

DFILE = os.path.join("{workingdir}", "_HerrMet_D{niter:03d}.npy")
DFILES = os.path.join("{workingdir}", "_HerrMet_D[0-9][0-9][0-9].npy")

FDFILE = os.path.join("{workingdir}", "_HerrMet_G{niter:03d}.npz")
FDFILES = os.path.join("{workingdir}", "_HerrMet_G[0-9][0-9][0-9].npz")

DOBSFILE = os.path.join("{workingdir}", "_HerrMet_Dobs.npy")
CDFILE = os.path.join("{workingdir}", "_HerrMet_CD.npz")
CDINVFILE = os.path.join("{workingdir}", "_HerrMet_CDinv.npz")
CMFILE = os.path.join("{workingdir}", "_HerrMet_CM.npz")
CMINVFILE = os.path.join("{workingdir}", "_HerrMet_CMinv.npz")


CLEAR_COMMAND = 'trash ' + " ".join(
    [NODEFILELOCAL.format(workingdir=WORKINGDIR),
     MPRIORFILE.format(workingdir=WORKINGDIR),
     DOBSFILE.format(workingdir=WORKINGDIR),
     CDFILE.format(workingdir=WORKINGDIR),
     CDINVFILE.format(workingdir=WORKINGDIR),
     CMFILE.format(workingdir=WORKINGDIR),
     CMINVFILE.format(workingdir=WORKINGDIR),
     MFILES.format(workingdir=WORKINGDIR),
     DFILES.format(workingdir=WORKINGDIR),
     FDFILES.format(workingdir=WORKINGDIR)])


def mfile2niter(mfile):
    return int(mfile.split('_M')[-1].split('.npy')[0])


def lastiter():
    path = MFILES.format(workingdir=WORKINGDIR)
    mfiles = glob.glob(path)
    if not len(mfiles):
        raise ValueError('no file responding to {}'.format(path))
    niters = [mfile2niter(mfile) for mfile in mfiles]
    return max(niters)


# ------------------------------ defaults
default_option = None

# ------------------------------ autorized_keys
authorized_keys = ["-option",
                   "-init",
                   "-prior",
                   "-data",
                   "-fd",
                   "-upd",
                   "-h", "-help"]

# ------------------------------ help messages
short_help = "--optimize   optimize the output of a point-wise inversion in 3D"

long_help = """\
--optimize   
    -init  s         clear temporary files from previous run, 
                     copy the input nodefile to {nodefilelocal} 
    -prior f f       get the prior model and model covariance matrix, 
                     save it to {mpriorfile}, {cmfile}
                     need the horizontal smoothing distance and vertical smoothing distance in km
    -data            get the data array and data covariance matrix, 
                     save it to {dobsfile}, {cdfile}
    -fd              compute the frechet derivatives for the first model
                     save it to {fdfile}                
    -h, -help        display the help message for this plugin 
""".format(workingdir=WORKINGDIR,
           nodefilelocal=NODEFILELOCAL.format(workingdir=WORKINGDIR),
           mpriorfile=MPRIORFILE.format(workingdir=WORKINGDIR),
           cmfile=CMFILE.format(workingdir=WORKINGDIR),
           dobsfile=DOBSFILE.format(workingdir=WORKINGDIR),
           cdfile=CDFILE.format(workingdir=WORKINGDIR),
           fdfile=FDFILES.format(workingdir=WORKINGDIR))

# ------------------------------ example usage
example = """\
## OPTIMIZE

HerrMet --optimize \\
    -init /home/max/prog/git/srfpython/tutorials/02_cube_inversion_example/nodes.txt \\
    -prior 50. 2.  \\
    -data \\
    -fd

HerrMet --optimize -upd
HerrMet --optimize -upd
HerrMet --optimize -upd    
    
"""


def haversine(loni, lati, lonj, latj, R=6371.):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    source https://stackoverflow.com/questions/15736995/how-can-i-quickly-estimate-the-distance-between-two-latitude-longitude-points
    lone, late : coordinates of the first point, in degrees
    lons, lats : coordinates of the second point, in degrees
    :return distance in km

    consistent with Ll2DGAB
    """
    # convert decimal degrees to radians
    q = np.pi / 180.

    # haversine formula
    dlon = (loni - lonj)
    dlat = (lati - latj)
    a = np.sin(q * dlat / 2.) ** 2. + np.cos(q * latj) * np.cos(q * lati) * np.sin(q * dlon / 2.) ** 2.
    c = 2. * np.arcsin(np.sqrt(a))
    km = R * c
    return km


def check_keys(argv):
    for k in argv.keys():
        if k in ['main', "_keyorder"]:
            continue  # private keys

        if k not in authorized_keys:
            message = 'option {} is not recognized'.format(k)

            warnings.warn(message)  # disable checks for development, please reactivate exception afterwards
            #raise Exception(message)


def get_superparameterizer(nodefile):
    """
    construct parameterizers needed to define the theory function in each node
    :param nodefile:
    :return:
    """
    # write the header of the parameterizer for each node
    parameter_string_header = """
    #met NLAYER = {}
    #met TYPE = 'mZVSVPvsRHvp'
    #met VPvs = 'lambda VS: 0.9409 + 2.0947 * VS - 0.8206 * VS ** 2.0 + 0.2683 * VS ** 3.0 - 0.0251 * VS ** 4.0'
    #met RHvp = 'lambda VP: 1.6612 * VP - 0.4721 * VP ** 2.0 + 0.0671 * VP ** 3.0 - 0.0043 * VP ** 4.0 + 0.000106 * VP ** 5.0'
    #fld KEY     VINF          VSUP
    #unt []      []            []
    #fmt %s      %f            %f
    """.format(len(nodefile.ztop)).replace('    #', '#')

    for i in range(1, len(nodefile.ztop)):
        # force VINF=VSUP => means lock the depth of the interfaces in the theory operator
        parameter_string_header += "-Z{} {} {}\n".format(i, -nodefile.ztop[i], -nodefile.ztop[i])  # add locked depth interfaces

    parameterizer_strings = []
    for nnode, (node, lon, lat, targetfile, paramfile, medianfile, p16file, p84file) in enumerate(nodefile):
        dm = depthmodel_from_mod96(medianfile)
        vs = dm.vs.interp(nodefile.zmid)

        parameter_string = parameter_string_header

        for i in range(len(nodefile.ztop)):
            # SET VINF < VS extracted from pointwise inv < VSUP
            # such as parameterizer.MMEAN corresponds to the extracted vs
            parameter_string += "VS{} {} {}\n".format(i, vs[i]-0.01, vs[i]+0.01)
        parameterizer_strings.append(parameter_string)

    parameterizer_strings = np.asarray(parameterizer_strings, str)
    parameterizers = [load_paramfile(parameter_string, verbose=False)[0]
                      for parameter_string in parameterizer_strings]
    return SuperParameterizer(parameterizers)


def get_superdatacoder(nodefile):
    """

    :param nodefile:
    :return:
    """
    datacoder_strings = []
    for nnode, (node, lon, lat, targetfile, paramfile, medianfile, p16file, p84file) in enumerate(nodefile):

        with open(targetfile, 'r') as fid:
            datacoder_string = "".join(fid.readlines())
            datacoder_strings.append(datacoder_string)

    datacoder_strings = np.asarray(datacoder_strings, str)
    datacoders = [makedatacoder(datacoder_string, which=Datacoder_log) for datacoder_string in datacoder_strings]
    return SuperDatacoder(datacoders)


def get_CM(nodefile, horizontal_smoothing_distance, vertical_smoothing_distance):
    ztop = nodefile.ztop
    zmid = nodefile.zmid

    CM_row_ind = np.array([], int)
    CM_col_ind = np.array([], int)
    CM_data = np.array([], float)
    for nnode in range(len(nodefile)):
        # find the posterior unc at each depth on vs from the pointwise inversion

        vs84_n = depthmodel_from_mod96(nodefile.p84files[nnode]).vs
        vs16_n = depthmodel_from_mod96(nodefile.p16files[nnode]).vs
        vs84_n = depthmodel1D(ztop, vs84_n.interp(zmid))
        vs16_n = depthmodel1D(ztop, vs16_n.interp(zmid))
        vs_unc_n = 0.5 * (vs84_n.values - vs16_n.values)

        # determine the vertical correlation coeff in cell n
        if vertical_smoothing_distance > 0:
            rhonn = np.exp(-0.5 * ((ztop - ztop[:, np.newaxis]) / vertical_smoothing_distance) ** 2.)
            covnn = vs_unc_n[:, np.newaxis] * vs_unc_n * rhonn

            row_ind, col_ind = np.meshgrid(range(len(ztop)), range(len(ztop)))
            CM_row_ind = np.hstack((CM_row_ind, nnode * len(ztop) + row_ind.flat[:]))
            CM_col_ind = np.hstack((CM_col_ind, nnode * len(ztop) + col_ind.flat[:]))
            CM_data = np.hstack((CM_data, covnn.flat[:]))
        else:
            covnn_diags = vs_unc_n ** 2.0
            CM_row_ind = np.hstack((CM_row_ind, nnode * len(ztop) + np.arange(len(ztop))))
            CM_col_ind = np.hstack((CM_col_ind, nnode * len(ztop) + np.arange(len(ztop))))
            CM_data = np.hstack((CM_data, covnn_diags.flat[:]))

        for mnode in range(nnode + 1, len(nodefile)):
            lonn = nodefile.lons[nnode]
            latn = nodefile.lats[nnode]
            lonm = nodefile.lons[mnode]
            latm = nodefile.lats[mnode]

            vs84_m = depthmodel_from_mod96(nodefile.p84files[mnode]).vs
            vs16_m = depthmodel_from_mod96(nodefile.p16files[mnode]).vs
            vs84_m = depthmodel1D(ztop, vs84_m.interp(zmid))
            vs16_m = depthmodel1D(ztop, vs16_m.interp(zmid))
            vs_unc_m = 0.5 * (vs84_m.values - vs16_m.values)

            dnm = haversine(lonn, latn, lonm, latm)
            rhonm_diags = np.exp(-0.5 * (dnm / horizontal_smoothing_distance) ** 2.)
            covnm_diags = vs_unc_n * vs_unc_m * rhonm_diags

            CM_row_ind = np.hstack((CM_row_ind, nnode * len(ztop) + np.arange(len(ztop))))
            CM_col_ind = np.hstack((CM_col_ind, mnode * len(ztop) + np.arange(len(ztop))))
            CM_data = np.hstack((CM_data, covnm_diags.flat[:]))

            CM_row_ind = np.hstack((CM_row_ind, mnode * len(ztop) + np.arange(len(ztop))))
            CM_col_ind = np.hstack((CM_col_ind, nnode * len(ztop) + np.arange(len(ztop))))
            CM_data = np.hstack((CM_data, covnm_diags.flat[:]))

    CM = csc_matrix((CM_data, (CM_row_ind, CM_col_ind)), shape=(len(nodefile) * len(ztop), len(nodefile) * len(ztop)))
    return CM


def get_CMinv(cmfile, cminvfile, verbose):
    try:
        CMinv = load_sparse_npz(cminvfile)
    except IOError:
        CM = load_sparse_npz(cmfile)
        CMinv = sparse_inv(CM)
        del CM
        save_matrix(cminvfile, CMinv, verbose)
    return CMinv


def save_matrix(filename, matrix, verbose):
    if verbose:
        print('saving {}, type:{}, shape:{}, dtype:{}'.format(
            filename, matrix.__class__.__name__, matrix.shape, matrix.dtype))

    if isinstance(matrix, np.ndarray):
        if not filename.endswith('.npy'):
            raise ValueError('{} does not end with .npy'.format(filename))
        np.save(filename, matrix, allow_pickle=False)

    elif isinstance(matrix, spmatrix):
        if not filename.endswith('.npz'):
            raise ValueError('{} does not end with .npz'.format(filename))
        save_sparse_npz(filename, matrix)

    else:
        raise NotImplementedError

    assert os.path.isfile(filename)


def get_propblem_size():
    n = len(np.load(DOBSFILE.format(workingdir=WORKINGDIR)))
    m = len(np.load(MPRIORFILE.format(workingdir=WORKINGDIR)))
    return n, m


class SuperParameterizer(object):
    def __init__(self, parameterizers):
        self.parameterizers = parameterizers

    def get_Mprior(self):
        Mprior = []
        for parameterizer in self.parameterizers:
            mprior = parameterizer.MMEAN
            Mprior = np.hstack((Mprior, mprior))
        return np.asarray(Mprior, float)

    def split(self, Model):
        """
        split a unique Model array into subarrays, one per node
        :param Model:
        :return:
        """
        j_current = 0
        models = []
        for node_number, parameterizer in enumerate(self.parameterizers):
            nparameters = len(parameterizer.MMEAN)
            model = Model[j_current : j_current + nparameters]
            j_current += nparameters
            models.append(model)

        assert j_current == len(Model)
        return models

    def inv(self, Model):
        models = []
        for nnode, model in enumerate(self.split(Model)):
            models.append(self.parameterizers[nnode].inv(model))
        return models

    def inv_to_depthmodels(self, Model):
        dms = []
        for nnode, model in enumerate(self.split(Model)):
            dms.append(self.parameterizers[nnode].inv_to_depthmodel(model))
        return dms

    def show_models(self, ax, Model, *args, **kwargs):
        for nnode, dm in enumerate(self.inv_to_depthmodels(Model)):
            dm.vs.values += nnode  # offset
            dm.vs.show(ax, *args, **kwargs)


class SuperDatacoder(object):
    def __init__(self, datacoders):
        self.datacoders = datacoders

    def get_Dobs(self):
        Dobs = []
        for datacoder in self.datacoders:
            dobs_current, CDinv_diag_current = datacoder.target()
            Dobs = np.concatenate((Dobs, dobs_current))
        return np.asarray(Dobs, float)

    def get_CD(self):
        CDinv_diag = []
        for datacoder in self.datacoders:
            dobs_node, CDinv_diag_node = datacoder.target()
            CDinv_diag.append(CDinv_diag_node)

        CDinv_diag = np.concatenate(CDinv_diag)

        CDinv = diags(CDinv_diag, offsets=0, format="csr", dtype=float)
        CD = diags(CDinv_diag ** -1.0, offsets=0, format="csr", dtype=float)
        return CD, CDinv

    def split(self, Data):
        i_current = 0
        datas = []
        for node_number, datacoder in enumerate(self.datacoders):
            ndatapoints = len(datacoder.values)
            data = Data[i_current : i_current + ndatapoints]
            i_current += ndatapoints
            datas.append(data)

        assert i_current == len(Data)
        return datas

    def inv(self, Data):
        datas = []
        for nnode, (datacoder, data) in enumerate(zip(self.datacoders, self.split(Data))):
            datas.append(datacoder.inv(data))

    def inv_to_laws(self, Data):
        Laws = []
        for nnode, (datacoder, data) in enumerate(zip(self.datacoders, self.split(Data))):
            Laws.append(datacoder.inv_to_laws(data))
        return Laws

    def show_datas(self, ax, Data, *args, **kwargs):
        for nnode, laws in enumerate(self.inv_to_laws(Data)):
            for law in laws:
                law.value += nnode
                law.show(ax, *args, **kwargs)


class SuperTheory(object):
    def __init__(self, superparameterizer, superdatacoder):
        self.superparameterizer = superparameterizer
        self.superdatacoder = superdatacoder
        self.theorys = []
        for parameterizer, datacoder in zip(superparameterizer.parameterizers,
                                            superdatacoder.datacoders):
            theory = Theory(parameterizer=parameterizer, datacoder=datacoder)
            self.theorys.append(theory)
        # number of data points for each node
        self.nds = np.asarray([len(theory.datacoder.values) for theory in self.theorys], int)

        # number of parameters for each node
        self.nps = np.asarray([len(theory.parameterizer.MMEAN) for theory in self.theorys], int)

    @property
    def shape(self):
        return self.nds.sum(), self.nps.sum()

    def __call__(self, Model):
        models = self.superparameterizer.split(Model)
        Data = []
        for theory, model in zip(self.theorys, models):
            data = theory(model)
            Data.append(data)
        return np.concatenate(Data)

    def get_FD(self, Model, mapkwargs):

        # prepare the frechet derivatives matrix
        G_shape = self.shape
        G_rows = np.hstack((0, self.nds.cumsum()[:-1]))  # index of the first row of G for each node
        G_cols = np.hstack((0, self.nps.cumsum()[:-1]))  # index of the first col of G for each node

        # compute the frechet derivatives for each node in parallel
        def node_generator():
            if Model is None:
                # means compute frechet derivatives relative to the prior model stored in p.MMEAN
                models = [theory.parameterizer.MMEAN for theory in self.theorys]
            else:
                # split the Model into sub models for each node
                models = self.superparameterizer.split(Model)

            for node_number, (theory, node_model) in enumerate(zip(self.theorys, models)):
                G_first_row, G_first_col = G_rows[node_number], G_cols[node_number]
                yield Job(node_number, theory, node_model, G_first_row, G_first_col)

        def node_calculator(node_number, node_theory, node_model, G_first_row, G_first_col):
            fd = node_theory.frechet_derivatives(m=node_model, gm=None)

            # compute the indices in fd
            fd_col_ind, fd_row_ind = np.meshgrid(range(fd.shape[1]), range(fd.shape[0]))

            # store indices and values to fill the sparse matrix G
            G_rows = G_first_row + fd_row_ind.flat[:]  # indices in the full matrix
            G_cols = G_first_col + fd_col_ind.flat[:]
            G_datas = fd.flat[:]

            return node_number, G_rows, G_cols, G_datas

        G_row_ind = np.array([], int)
        G_col_ind = np.array([], int)
        G_fd_data = np.array([], float)
        with MapSync(node_calculator, node_generator(), **mapkwargs) as ma:
            for jobid, (node_number, G_rows, G_cols, G_datas), _, _ in ma:
                G_row_ind = np.concatenate((G_row_ind, G_rows))
                G_col_ind = np.concatenate((G_col_ind, G_cols))
                G_fd_data = np.concatenate((G_fd_data, G_datas))

        G = csr_matrix((G_fd_data, (G_row_ind, G_col_ind)), shape=G_shape)
        return G


def optimize(argv, verbose, mapkwargs):
    # import matplotlib.pyplot as plt
    # import numpy as np
    #
    # mprior = np.load('_HerrMet_Mprior.npy')
    # m000 = np.load('_HerrMet_M000.npy')
    # m001 = np.load('_HerrMet_M001.npy')
    # m002 = np.load('_HerrMet_M002.npy')
    # plt.plot(m000, "k")
    # plt.plot(m001, "g")
    # plt.plot(m002, "b")
    # plt.show()
    #
    # dobs = np.load('_HerrMet_Dobs.npy')
    # plt.plot(dobs, "k")
    # d000 = np.load('_HerrMet_D000.npy')
    # plt.plot(d000, "g")
    # plt.plot(d000 - dobs, "g")
    # d001 = np.load('_HerrMet_D001.npy')
    # plt.plot(d001, "b")
    # plt.plot(d001 -dobs, "b")
    # d002 = np.load('_HerrMet_D002.npy')
    # plt.plot(d002, "r")
    # plt.plot(d002-dobs, "r")
    # plt.show()
    #
    #
    # sys.exit()
    if '-h' in argv.keys() or "-help" in argv.keys():
        print(long_help)
        return

    check_keys(argv)

    # === inititate the working directory, clean up previous files, copy node file
    if "-init" in argv.keys():
        nodefilename_in = argv['-init'][0]

        # ========= clear temporary files
        print(CLEAR_COMMAND)
        os.system(CLEAR_COMMAND)

        # ========= get the node file and copy it in the working directory
        if not os.path.isfile(nodefilename_in):
            raise IOError(nodefilename_in)

        nodefilename_out = NODEFILELOCAL.format(workingdir=WORKINGDIR)

        nodefile = NodeFile(nodefilename_in)
        nodefile.fill_extraction_files()
        nodefile.write(nodefilename_out)
        assert os.path.isfile(nodefilename_out)

    # =================================
    # (re)read the local node file
    nodefile = NodeFile(NODEFILELOCAL.format(workingdir=WORKINGDIR))
    nodefile.fill_extraction_files()

    # Compute the initialization matrixes
    superparameterizer = get_superparameterizer(nodefile)
    superdatacoder = get_superdatacoder(nodefile)
    supertheory = SuperTheory(
        superparameterizer=superparameterizer,
        superdatacoder=superdatacoder)
    n_data_points, n_parameters = supertheory.shape

    if "-show" in argv.keys():
        fig_model = plt.figure()
        ax_model = fig_model.add_subplot(111)
        modelfiles = np.sort(glob.glob(MFILES.format(workingdir=WORKINGDIR)))
        for modelfile, color in zip([modelfiles[0], modelfiles[-1]], "kr"):
            Model = np.load(modelfile)
            superparameterizer.show_models(ax_model, Model, color=color)

        fig_data = plt.figure()
        ax_data = fig_data.add_subplot(111)
        dobsfile = DOBSFILE.format(workingdir=WORKINGDIR)
        Dobs = np.load(dobsfile)
        superdatacoder.show_datas(ax_data, Dobs, color="k", showdvalue=True)

        datafiles = np.sort(glob.glob(DFILES.format(workingdir=WORKINGDIR)))
        for datafile, color in zip([datafiles[0], datafiles[-1]], 'br'):
            Data = np.load(datafile)
            superdatacoder.show_datas(ax_data, Data, color=color, showdvalue=False)

        plt.ion()
        plt.show()
        input('pause')
        exit()

    if "-prior" in argv.keys():
        hsd, vsd = argv["-prior"]
        # ========= set the prior model and covariance matrix
        mpriorfilename = MPRIORFILE.format(workingdir=WORKINGDIR)
        cmfilename = CMFILE.format(workingdir=WORKINGDIR)
        Mprior = superparameterizer.get_Mprior()
        CM = get_CM(nodefile,
                    horizontal_smoothing_distance=hsd,
                    vertical_smoothing_distance=vsd)
        save_matrix(mpriorfilename, Mprior, verbose)
        save_matrix(cmfilename, CM, verbose)

        M0 = Mprior
        m0filename = MFILE.format(workingdir=WORKINGDIR, niter=0)
        save_matrix(m0filename, M0, verbose)

        d0filename = DFILE.format(workingdir=WORKINGDIR, niter=0)
        D0 = supertheory(M0)
        save_matrix(d0filename, D0, verbose)

    if "-data" in argv.keys():
        # ========= set the observed data and data covariance matrix
        Dobs = superdatacoder.get_Dobs()
        CD, CDinv = superdatacoder.get_CD()
        dobsfilename = DOBSFILE.format(workingdir=WORKINGDIR)
        cdfilename = CDFILE.format(workingdir=WORKINGDIR)
        cdinvfilename = CDINVFILE.format(workingdir=WORKINGDIR)

        save_matrix(dobsfilename, Dobs, verbose)
        save_matrix(cdfilename, CD, verbose)
        save_matrix(cdinvfilename, CDinv, verbose)

    if "-fd" in argv.keys():
        # ========= compute the frechet derivatives
        niter = lastiter()
        mfilename = MFILE.format(workingdir=WORKINGDIR, niter=niter)
        fdfilename = FDFILE.format(workingdir=WORKINGDIR, niter=niter)
        if not os.path.isfile(fdfilename):
            Mi = np.load(mfilename)
            Gi = supertheory.get_FD(Mi, mapkwargs)
            save_matrix(fdfilename, Gi, verbose)

        elif verbose:
            print('{} exists already'.format(fdfilename))

    if "-upd" in argv.keys():

        niter = lastiter()
        m0file = MFILE.format(workingdir=WORKINGDIR, niter=0)
        mfilename = MFILE.format(workingdir=WORKINGDIR, niter=niter)
        dfile = DFILE.format(workingdir=WORKINGDIR, niter=niter)
        fdfilename = FDFILE.format(workingdir=WORKINGDIR, niter=niter)
        cmfile = CMFILE.format(workingdir=WORKINGDIR)
        cminvfile = CMINVFILE.format(workingdir=WORKINGDIR)
        cdfile = CDFILE.format(workingdir=WORKINGDIR)
        cdinvfile = CDINVFILE.format(workingdir=WORKINGDIR)
        dobsfile = DOBSFILE.format(workingdir=WORKINGDIR)

        M0 = np.load(m0file)
        Mi = np.load(mfilename)
        Di = np.load(dfile)  # g(Mi)
        Gi = load_sparse_npz(fdfilename)
        Dobs = np.load(dobsfile)

        if n_data_points <= n_parameters:
            # over determined problem
            CM = load_sparse_npz(cmfile)
            CD = load_sparse_npz(cdfile)

            CMGiT = CM * Gi.T
            Siinv = sparse_inv(CD + Gi * CMGiT)
            Hi = np.dot(CMGiT, Siinv)

        else:
            CDinv = load_sparse_npz(cdinvfile)
            CMinv = get_CMinv(cmfile, cminvfile, verbose)

            GiTCDinv = Gi.T * CDinv
            Siinv = sparse_inv(GiTCDinv * Gi + CMinv)
            Hi = Siinv * GiTCDinv

        Xi = Dobs - Di + Gi * (Mi - M0)
        Minew = M0 + Hi * Xi

        Dinew = supertheory(Minew)
        Ginew = supertheory.get_FD(Model=Minew, mapkwargs=mapkwargs)

        mfilename_new = MFILE.format(workingdir=WORKINGDIR, niter=niter + 1)
        dfilename_new = DFILE.format(workingdir=WORKINGDIR, niter=niter + 1)
        fdfilename_new = FDFILE.format(workingdir=WORKINGDIR, niter=niter + 1)

        save_matrix(mfilename_new, Minew, verbose)
        save_matrix(dfilename_new, Dinew, verbose)
        save_matrix(fdfilename_new, Ginew, verbose)



