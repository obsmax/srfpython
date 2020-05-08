from __future__ import print_function
from builtins import input
import sys, glob, os
import warnings
import numpy as np
from scipy.sparse import spmatrix, diags, issparse, csc_matrix, \
    save_npz as save_sparse_npz, load_npz as load_sparse_npz
from scipy.sparse.linalg import spsolve, splu  #inv as sparse_inv NO!!!
import matplotlib.pyplot as plt
from srfpython.HerrMet.nodefile import NodeFile
from srfpython.depthdisp.depthmodels import depthmodel_from_mod96, depthmodel1D, depthmodel_from_arrays
from srfpython.HerrMet.paramfile import load_paramfile
from srfpython.HerrMet.datacoders import makedatacoder, Datacoder_log
from srfpython.HerrMet.theory import Theory
from srfpython.standalone.multipro8 import Job, MapSync
from srfpython.Herrmann.Herrmann import CPiSDomainError
from srfpython.HerrMet.files import ROOTKEY

WORKINGDIR = "."
NODEFILELOCAL = os.path.join(WORKINGDIR, ROOTKEY + "nodes.txt")
MPRIORFILE = os.path.join(WORKINGDIR, ROOTKEY + "Mprior.npy")
MFILE = os.path.join(WORKINGDIR, ROOTKEY + "M{niter:03d}.npy")
MFILES = os.path.join(WORKINGDIR, ROOTKEY + "M[0-9][0-9][0-9].npy")
DFILE = os.path.join(WORKINGDIR, ROOTKEY + "D{niter:03d}.npy")
DFILES = os.path.join(WORKINGDIR, ROOTKEY + "D[0-9][0-9][0-9].npy")
FDFILE = os.path.join(WORKINGDIR, ROOTKEY + "G{niter:03d}.npz")
FDFILES = os.path.join(WORKINGDIR, ROOTKEY + "G[0-9][0-9][0-9].npz")
DOBSFILE = os.path.join(WORKINGDIR, ROOTKEY + "Dobs.npy")
CDFILE = os.path.join(WORKINGDIR, ROOTKEY + "CD.npz")
CDINVFILE = os.path.join(WORKINGDIR, ROOTKEY + "CDinv.npz")
CMFILE = os.path.join(WORKINGDIR, ROOTKEY + "CM.npz")
CMINVFILE = os.path.join(WORKINGDIR, ROOTKEY + "CMinv.npz")


CLEAR_COMMAND = 'trash ' + " ".join(
    [NODEFILELOCAL,
     MPRIORFILE,
     DOBSFILE,
     CDFILE,
     CDINVFILE,
     CMFILE,
     MFILES,
     DFILES,
     FDFILES])


def mfile2niter(mfile):
    return int(mfile.split('_M')[-1].split('.npy')[0])


def lastiter():
    path = MFILES
    mfiles = glob.glob(path)
    if not len(mfiles):
        raise ValueError('no file responding to {}'.format(path))
    niters = [mfile2niter(mfile) for mfile in mfiles]
    return max(niters)


# ------------------------------ defaults
default_option = None

# ------------------------------ autorized_keys
authorized_keys = ["-option", "-init", "-restart", "-prior",
                   "-data", "-fd", "-upd", "-show",
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
    -upd   i         update the model, compute the frechet derivatives for the next iteration
                     repeat n times, default 1                                 
    -h, -help        display the help message for this plugin 
""".format(workingdir=WORKINGDIR,
           nodefilelocal=NODEFILELOCAL,
           mpriorfile=MPRIORFILE,
           cmfile=CMFILE,
           dobsfile=DOBSFILE,
           cdfile=CDFILE,
           fdfile=FDFILES)

# ------------------------------ example usage
example = """\
## OPTIMIZE

# initialize the optimization (clear temp files)
HerrMet --optimize \\
    -init /home/max/prog/git/srfpython/tutorials/02_cube_inversion_example/nodes.txt 
    
# set the data and prior information    
HerrMet --optimize \\
    -prior 80. 0.  \\
    -data 
    
# the prior can be changed afterwards, 
# it requires to remove the iteration files    
HerrMet --optimize \\
    -restart \\
    -prior 100. 0.

# set the sensituvuty kernels for first iteration
HerrMet --optimize -fd 
    
# run 3 iterations
HerrMet --optimize -upd 3 

# display results
HerrMet --optimize -show    
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

            # warnings.warn(message)  # disable checks for development, please reactivate exception afterwards
            raise Exception(message)


class NodeFileLocal(NodeFile):

    def get_superparameterizer(self):
        """
        construct parameterizers needed to define the theory function in each node
        :param self:
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
        """.format(len(self.ztop)).replace('    #', '#')

        for i in range(1, len(self.ztop)):
            # force VINF=VSUP => means lock the depth of the interfaces in the theory operator
            parameter_string_header += "-Z{} {} {}\n".format(i, -self.ztop[i], -self.ztop[i])  # add locked depth interfaces

        parameterizer_strings = []
        for nnode, (node, lon, lat, targetfile, paramfile, medianfile, p16file, p84file) in enumerate(self):
            dm = depthmodel_from_mod96(medianfile)
            vs = dm.vs.interp(self.zmid)

            parameter_string = parameter_string_header

            for i in range(len(self.ztop)):
                # SET VINF < VS extracted from pointwise inv < VSUP
                # such as parameterizer.MMEAN corresponds to the extracted vs
                parameter_string += "VS{} {} {}\n".format(i, vs[i]-0.01, vs[i]+0.01)
            parameterizer_strings.append(parameter_string)

        parameterizer_strings = np.asarray(parameterizer_strings, str)
        parameterizers = [load_paramfile(parameter_string, verbose=False)[0]
                          for parameter_string in parameterizer_strings]
        return SuperParameterizer(parameterizers)

    def get_superdatacoder(self):
        """

        :param self:
        :return:
        """
        datacoder_strings = []
        for nnode, (node, lon, lat, targetfile, paramfile, medianfile, p16file, p84file) in enumerate(self):

            with open(targetfile, 'r') as fid:
                datacoder_string = "".join(fid.readlines())
                datacoder_strings.append(datacoder_string)

        datacoder_strings = np.asarray(datacoder_strings, str)
        datacoders = [makedatacoder(datacoder_string, which=Datacoder_log) for datacoder_string in datacoder_strings]
        return SuperDatacoder(datacoders)

    def get_CM(self, horizontal_smoothing_distance, vertical_smoothing_distance):
        ztop = self.ztop
        zmid = self.zmid
        lock_half_space = True
        lock_half_space_sigma = 0.001

        CM_row_ind = np.array([], int)
        CM_col_ind = np.array([], int)
        CM_data = np.array([], float)

        for nnode in range(len(self)):
            # find the posterior unc at each depth on vs from the pointwise inversion

            vs84_n = depthmodel_from_mod96(self.p84files[nnode]).vs
            vs16_n = depthmodel_from_mod96(self.p16files[nnode]).vs
            vs84_n = depthmodel1D(ztop, vs84_n.interp(zmid))
            vs16_n = depthmodel1D(ztop, vs16_n.interp(zmid))
            vs_unc_n = 0.5 * (vs84_n.values - vs16_n.values)
            # vs_unc_n = np.ones(len(ztop))  ######################################" TEST
            if lock_half_space:
                vs_unc_n[-1] = lock_half_space_sigma

            # determine the vertical correlation coeff in cell n
            if vertical_smoothing_distance > 0:
                raise Exception('seems incorrect (check the sign of (M - Mprior)T * CMinv * (M - Mprior)')
                rhonn = np.exp(-0.5 * ((ztop - ztop[:, np.newaxis]) / vertical_smoothing_distance) ** 2.)
                covnn = vs_unc_n[:, np.newaxis] * vs_unc_n * rhonn
                col_ind, row_ind = np.meshgrid(range(len(ztop)),
                                               range(len(ztop)))
                CM_row_ind = np.hstack((CM_row_ind, nnode * len(ztop) + row_ind.flat[:]))
                CM_col_ind = np.hstack((CM_col_ind, nnode * len(ztop) + col_ind.flat[:]))
                CM_data = np.hstack((CM_data, covnn.flat[:]))

            else:
                covnn_diags = vs_unc_n ** 2.0
                CM_row_ind = np.hstack((CM_row_ind, nnode * len(ztop) + np.arange(len(ztop))))
                CM_col_ind = np.hstack((CM_col_ind, nnode * len(ztop) + np.arange(len(ztop))))
                CM_data = np.hstack((CM_data, covnn_diags.flat[:]))

            if horizontal_smoothing_distance > 0:
                for mnode in range(nnode + 1, len(self)):
                    lonn = self.lons[nnode]
                    latn = self.lats[nnode]
                    lonm = self.lons[mnode]
                    latm = self.lats[mnode]

                    vs84_m = depthmodel_from_mod96(self.p84files[mnode]).vs
                    vs16_m = depthmodel_from_mod96(self.p16files[mnode]).vs
                    vs84_m = depthmodel1D(ztop, vs84_m.interp(zmid))
                    vs16_m = depthmodel1D(ztop, vs16_m.interp(zmid))
                    vs_unc_m = 0.5 * (vs84_m.values - vs16_m.values)
                    # vs_unc_m = np.ones(len(ztop))  ######################################" TEST

                    if lock_half_space:
                        vs_unc_m[-1] = lock_half_space_sigma

                    dnm = haversine(loni=lonn, lati=latn, lonj=lonm, latj=latm)
                    rhonm = np.exp(-0.5 * (dnm / horizontal_smoothing_distance) ** 2.)
                    if True:
                        covnm = vs_unc_n * vs_unc_m * rhonm

                        CM_row_ind = np.hstack((CM_row_ind, nnode * len(ztop) + np.arange(len(ztop))))
                        CM_col_ind = np.hstack((CM_col_ind, mnode * len(ztop) + np.arange(len(ztop))))
                        CM_data = np.hstack((CM_data, covnm.flat[:]))

                        CM_row_ind = np.hstack((CM_row_ind, mnode * len(ztop) + np.arange(len(ztop))))
                        CM_col_ind = np.hstack((CM_col_ind, nnode * len(ztop) + np.arange(len(ztop))))
                        CM_data = np.hstack((CM_data, covnm.flat[:]))

                    else:
                        covnm = vs_unc_n[:, np.newaxis] * vs_unc_m * rhonm

                        col_ind, row_ind = np.meshgrid(range(len(ztop)), range(len(ztop)))

                        CM_row_ind = np.hstack((CM_row_ind, nnode * len(ztop) + row_ind.flat[:]))
                        CM_col_ind = np.hstack((CM_col_ind, mnode * len(ztop) + col_ind.flat[:]))
                        CM_data = np.hstack((CM_data, covnm.flat[:]))

                        CM_row_ind = np.hstack((CM_row_ind, mnode * len(ztop) + col_ind.flat[:]))
                        CM_col_ind = np.hstack((CM_col_ind, nnode * len(ztop) + row_ind.flat[:]))
                        CM_data = np.hstack((CM_data, covnm.flat[:]))

        CM = csc_matrix((CM_data, (CM_row_ind, CM_col_ind)),
                        shape=(len(self) * len(ztop), len(self) * len(ztop)))

        check = (CM - CM.T)
        assert not len(check.data)
        # print(check.format, check.shape, check.data)

        return CM


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

    def show_models(self, ax, Model, Munc=None, offset=3.0, **kwargs):
        xticks = []
        xticklabels = []
        vsticks = np.array([1., 2., 3.])
        vsticklabels = ['1', '2', '3']
        if Munc is None:
            for nnode, dm in enumerate(self.inv_to_depthmodels(Model)):
                dm.vs.values += offset * nnode  # offset
                dm.vs.show(ax, **kwargs)
                xticks = np.concatenate((xticks, vsticks + offset * nnode))
                xticklabels = np.concatenate((xticklabels, vsticklabels))
        else:
            for nnode, (dmlow, dmhigh) in enumerate(zip(self.inv_to_depthmodels(Model-Munc),
                                                        self.inv_to_depthmodels(Model+Munc))):
                dmlow.vs.values += offset * nnode  # offset
                dmhigh.vs.values += offset * nnode  # offset
                dmlow.vs.fill_between(ax=ax, other=dmhigh.vs, **kwargs)
                xticks = np.concatenate((xticks, vsticks + offset * nnode))
                xticklabels = np.concatenate((xticklabels, vsticklabels))
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)


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

        CDinv = diags(CDinv_diag, offsets=0, format="csc", dtype=float)
        CD = diags(CDinv_diag ** -1.0, offsets=0, format="csc", dtype=float)
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

    def show_datas(self, ax, Data, showdvalue=False,
                   offset=3.,
                   pticks=[0.5, 1., 2.],
                   pticklabels=['0.5', '1', '2'],
                   **kwargs):
        xticks = []
        xticklabels = []

        for nnode, laws in enumerate(self.inv_to_laws(Data)):

            xticks = np.concatenate((xticks, np.log(pticks) + offset * nnode))
            xticklabels = np.concatenate((xticklabels, pticklabels))
            for law in laws:

                x = np.log(1. / law.freq) + offset * nnode
                y = law.value
                if showdvalue:
                    yinf = law.value * np.exp(-law.dvalue / law.value)
                    ysup = law.value * np.exp(+law.dvalue / law.value)
                    ax.fill_between(x, yinf, ysup, **kwargs)
                else:
                    ax.plot(x, y, **kwargs)
        ax.set_xscale('linear')
        ax.set_yscale('log')
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)


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

        G = csc_matrix((G_fd_data, (G_row_ind, G_col_ind)), shape=G_shape)
        return G


def sparsedet(M):
    lu = splu(M)
    diagL = lu.L.diagonal()
    diagU = lu.U.diagonal()
    diagL = diagL.astype(np.complex128)
    diagU = diagU.astype(np.complex128)
    logdet = np.log(diagL).sum() + np.log(diagU).sum()
    det = np.exp(logdet)  # usually underflows/overflows for large matrices
    return det.real


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


def Ainv_dot_b(A, b):
    """
    computes A^-1 * b without inverting A
    for numerical stability
    see Tarantola 2005, section 3.4.5 p80
    :param A:
    :param b:
    :return:
    """
    if issparse(A):
        return spsolve(A=A, b=b)

    elif isinstance(A, np.ndarray) or isinstance(A, np.matrix):
        return np.linalg.solve(a=A, b=b)

    else:
        raise TypeError


def optimize(argv, verbose, mapkwargs):

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

        nodefile = NodeFile(nodefilename_in)
        nodefile.fill_extraction_files()
        nodefile.write(NODEFILELOCAL)
        assert os.path.isfile(NODEFILELOCAL)

    if "-restart" in argv.keys():
        # ========= clear temporary files except for the initiation files
        for niter in range(1, lastiter()+1):
            mfile = MFILE.format(niter=niter)
            dfile = DFILE.format(niter=niter)
            fdfile = FDFILE.format(niter=niter)
            trashcmd = 'trash {} {} {}'.format(mfile, dfile, fdfile)
            print(trashcmd)
            os.system(trashcmd)

    # =================================
    # (re)read the local node file
    nodefile = NodeFileLocal(NODEFILELOCAL)
    nodefile.fill_extraction_files()

    # Compute the initialization matrixes
    superparameterizer = nodefile.get_superparameterizer()
    superdatacoder = nodefile.get_superdatacoder()
    supertheory = SuperTheory(
        superparameterizer=superparameterizer,
        superdatacoder=superdatacoder)
    n_data_points, n_parameters = supertheory.shape

    if "-prior" in argv.keys():
        hsd, vsd = argv["-prior"]
        # ========= set the prior model and covariance matrix

        Mprior = superparameterizer.get_Mprior()
        CM = nodefile.get_CM(
            horizontal_smoothing_distance=hsd,
            vertical_smoothing_distance=vsd)
        save_matrix(MPRIORFILE, Mprior, verbose)
        save_matrix(CMFILE, CM, verbose)

        M0 = Mprior
        m0filename = MFILE.format(niter=0)
        save_matrix(m0filename, M0, verbose)

        d0filename = DFILE.format(niter=0)
        D0 = supertheory(M0)
        save_matrix(d0filename, D0, verbose)

    if "-data" in argv.keys():
        # ========= set the observed data and data covariance matrix
        Dobs = superdatacoder.get_Dobs()
        CD, CDinv = superdatacoder.get_CD()

        save_matrix(DOBSFILE, Dobs, verbose)
        save_matrix(CDFILE, CD, verbose)
        save_matrix(CDINVFILE, CDinv, verbose)

    if "-fd" in argv.keys():
        # ========= compute the frechet derivatives for the last model found
        niter = lastiter()
        mfilename = MFILE.format(niter=niter)
        fdfilename = FDFILE.format(niter=niter)
        if not os.path.isfile(fdfilename):
            Mi = np.load(mfilename)
            Gi = supertheory.get_FD(Mi, mapkwargs)
            save_matrix(fdfilename, Gi, verbose)

        elif verbose:
            print('{} exists already'.format(fdfilename))

    if "-upd" in argv.keys():
        nupd = argv["-upd"][0] if len(argv["-upd"]) > 0 else 1

        # ================ load data
        Mprior = np.load(MPRIORFILE)
        M0 = np.load(MFILE.format(niter=0))
        Dobs = np.load(DOBSFILE)

        CM = load_sparse_npz(CMFILE)
        CD = load_sparse_npz(CDFILE)

        for _ in range(nupd):
            niter = lastiter()

            # ==== save the current model state
            mfilename = MFILE.format(niter=niter)
            dfile = DFILE.format(niter=niter)
            fdfilename = FDFILE.format(niter=niter)

            Mi = np.load(mfilename)
            Di = np.load(dfile)  # g(Mi)
            Gi = load_sparse_npz(fdfilename)

            # ==== update the current model
            has_nan = np.isnan(Dobs).any()
            Dinew, Minew = tv23_1(Dobs, Di, CD,
                       Mprior, M0, Mi, CM,
                       Gi, supertheory)
            if has_nan:
                assert np.isnan(Dobs).any()

            # ==== compute the frechet derivatives for the next iteration
            Ginew = supertheory.get_FD(Model=Minew, mapkwargs=mapkwargs)

            # ==== save the new model state
            mfilename_new = MFILE.format(niter=niter + 1)
            dfilename_new = DFILE.format(niter=niter + 1)
            fdfilename_new = FDFILE.format(niter=niter + 1)

            save_matrix(mfilename_new, Minew, verbose)
            save_matrix(dfilename_new, Dinew, verbose)
            save_matrix(fdfilename_new, Ginew, verbose)

    if "-show" in argv.keys():
        mpriorfile = MPRIORFILE
        modelfiles = np.sort(glob.glob(MFILES))
        datafiles = np.sort(glob.glob(DFILES))

        Dobs = np.load(DOBSFILE)
        Mprior = np.load(mpriorfile)
        CM = load_sparse_npz(CMFILE)
        CDinv = load_sparse_npz(CDINVFILE)
        Mpriorunc = np.asarray(CM[np.arange(n_parameters), np.arange(n_parameters)]).flat[:] ** 0.5

        # =================== display the data fit
        fig_costs = plt.figure()
        ax_costs1 = fig_costs.add_subplot(311)
        ax_costs2 = fig_costs.add_subplot(312, sharex=ax_costs1)
        ax_costs3 = fig_costs.add_subplot(313, sharex=ax_costs1)

        data_costs = []
        model_costs = []
        for niter in range(lastiter()+1):
            modelfile = MFILE.format(niter=niter)
            datafile = DFILE.format(niter=niter)
            Model = np.load(modelfile)
            Data = np.load(datafile)

            assert issparse(CDinv)
            assert issparse(CM)
            assert not issparse(Data - Dobs)
            assert not issparse(Model - Mprior)

            # nan issues
            # nan can occur in the predicted data if a data point is above the cut off period
            # let ignore them
            Inan = np.isnan(Data) | np.isnan(Dobs)
            Dobs[Inan] = Data[Inan] = 0.  # so that they do not act on the data misfit

            data_cost = np.dot((Data - Dobs), CDinv * (Data - Dobs))  # only because CDinv is known exactly
            model_cost = np.dot(Model - Mprior, Ainv_dot_b(A=CM, b=Model - Mprior))
            if verbose:
                print('iter {:03d}: chi2_data={:<16f} chi2_model={:<16f} chi2={:<16f}'.format(
                    niter, data_cost, model_cost, data_cost + model_cost))
            # model_cost = np.dot(np.dot((Model - Mprior), CMinv), (Model - Mprior))  # NO
            data_costs.append(data_cost)
            model_costs.append(model_cost)

        ax_costs1.plot(data_costs, 'ko')
        ax_costs2.plot(model_costs, 'bo')
        ax_costs3.plot(np.array(data_costs) + np.array(model_costs), 'bo')

        # =================== display the models side by side
        fig_model = plt.figure()
        ax_model = fig_model.add_subplot(111)

        superparameterizer.show_models(ax_model, Mprior - Mpriorunc, color="k", alpha=0.3)
        superparameterizer.show_models(ax_model, Mprior + Mpriorunc, color="k", alpha=0.3)
        superparameterizer.show_models(ax_model, Mprior, color="k", alpha=1.0)
        superparameterizer.show_models(ax_model, Mprior, Munc=Mpriorunc, color="k", alpha=0.3)

        for modelfile in modelfiles[1:-1]:
            Model = np.load(modelfile)
            superparameterizer.show_models(ax_model, Model, color="k", alpha=0.1)

        for modelfile, color in zip([modelfiles[0], modelfiles[-1]], "br"):
            Model = np.load(modelfile)
            superparameterizer.show_models(ax_model, Model, color=color, linewidth=2)

        # =================== display the datas side by side
        fig_data = plt.figure()
        ax_data = fig_data.add_subplot(111)

        superdatacoder.show_datas(ax_data, Dobs, color="k", showdvalue=True, alpha=0.3)
        superdatacoder.show_datas(ax_data, Dobs, color="k", showdvalue=False)

        for datafile, color in zip([datafiles[0], datafiles[-1]], 'br'):
            Data = np.load(datafile)
            superdatacoder.show_datas(ax_data, Data, color=color, showdvalue=False)

        plt.ion()
        plt.show()
        input('pause')
        exit()


def tv23_1(Dobs, Di, CD,
           Mprior, M0, Mi, CM,
           Gi, supertheory):

    # nan issues
    # nan can occur in the predicted data if a data point is above the cut off period
    # let ignore them
    Inan = np.isnan(Di) | np.isnan(Dobs)
    Dobs[Inan] = Di[Inan] = 0.  # TODO make sure this has not affected Dobs outside this function
    Xi = Dobs - Di + Gi * (Mi - Mprior)

    CMGiT = CM * Gi.T
    KiXi = CMGiT * Ainv_dot_b(A=CD + Gi * CMGiT, b=Xi)

    mu = 1.0
    Minew = Dinew = None
    while mu > 0.001:
        try:
            Minew = Mi + mu * (M0 + KiXi - Mi)
            Dinew = supertheory(Minew)
            break
        except CPiSDomainError as e:
            print('step was too large and leaded to an error of the forward problem, '
                  'try reducing the step')
            mu /= 2.0
            continue

    if Minew is None:
        raise Exception('fail')

    return Dinew, Minew
