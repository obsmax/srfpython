from __future__ import print_function
from builtins import input
import warnings
import glob, os
import numpy as np
from scipy.sparse import spmatrix, diags, issparse, csc_matrix, \
    save_npz as save_sparse_npz, load_npz as load_sparse_npz
from scipy.sparse.linalg import spsolve, splu  #inv as sparse_inv NO!!!
import matplotlib.pyplot as plt
from srfpython.HerrMet.nodefile import NodeFile
from srfpython.depthdisp.depthmodels import depthmodel_from_mod96, depthmodel1D
from srfpython.HerrMet.paramfile import load_paramfile
from srfpython.HerrMet.datacoders import makedatacoder, Datacoder_log
from srfpython.HerrMet.theory import Theory
from srfpython.standalone.multipro8 import Job, MapSync
from srfpython.Herrmann.Herrmann import CPiSDomainError
from srfpython.HerrMet.files import ROOTKEY
from srfpython.coordinates import haversine

# TODO replace shell commands using python
# TODO add convergence interruption

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

M96FILEOUT = os.path.join(WORKINGDIR, ROOTKEY + "{node}_optimized.mod96")
S96FILEOUT = os.path.join(WORKINGDIR, ROOTKEY + "{node}_optimized.surf96")
M96FILEOUTS = os.path.join(WORKINGDIR, ROOTKEY + "*_optimized.mod96")
S96FILEOUTS = os.path.join(WORKINGDIR, ROOTKEY + "*_optimized.surf96")

CLEAR_COMMAND = 'trash ' + " ".join(
    [NODEFILELOCAL,
     MPRIORFILE,
     DOBSFILE,
     CDFILE,
     CDINVFILE,
     CMFILE,
     MFILES,
     DFILES,
     FDFILES,
     M96FILEOUTS,
     S96FILEOUTS])


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
                   "-data", "-fd", "-upd", "-show", "-save",
                   "-h", "-help"]

# ------------------------------ help messages
short_help = "--optimize   optimize the output of a point-wise inversion in 3D"

long_help = """\
--optimize   
    -init  s         clear temporary files from previous run, 
                     copy the input nodefile to {nodefilelocal} 
    -prior fx6       get the prior model and model covariance matrix, 
                     save it to {mpriorfile}, {cmfile}
                     needs 6 arguments
                        - horizontal_smoothing_distance in km 
                        - vertical_smoothing_distance in km 
                        - trunc_horizontal_smoothing_distance in km 
                        - trunc_vertical_smoothing_distance in km 
                        - lock_half_space 1 or 0 
                        - scale_uncertainties (give a factor >= 0) 
                        - add_uncertainty in km/s (>= 0)
    -data  fx2       get the data array and data covariance matrix, 
                     save it to {dobsfile}, {cdfile}
                     needs 2 arguments
                        - scale_uncertainties (give a factor >= 0) 
                        - add_uncertainty in km/s (>= 0)
    -fd              compute the frechet derivatives for last iteration found
                     save it to {fdfile}
    -upd  [i [i]]    run iterative inversion
                     two optional arguments
                        - provide the number of iterations to run (default 1)
                        - update the frechet derivatives before inversion 1/0 (default 1)
                          you can also not update the fr. der. for several iterations
                          and update it manually using -fd                                 
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
    -prior 100. 2. 1000. 20. 1 1.0 0.5   \\
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
        #met VPvs = 'lambda VS: 0.9409+2.0947*VS-0.8206*VS**2+0.2683*VS**3-0.0251*VS**4'
        #met RHvp = 'lambda VP: 1.6612*VP-0.4721*VP**2+0.0671*VP**3-0.0043*VP**4+0.000106*VP**5'
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

    def get_CM(
            self, horizontal_smoothing_distance, vertical_smoothing_distance,
            trunc_horizontal_smoothing_distance=0.,
            trunc_vertical_smoothing_distance=0.,
            lock_half_space=True,
            lock_half_space_sigma=0.01,
            scale_uncertainties=1.,
            add_uncertainty=0.,
            norm="L1",
            visual_qc=False):

        ztop = self.ztop
        zmid = self.zmid

        if trunc_vertical_smoothing_distance is None or trunc_vertical_smoothing_distance <= 0.:
            # truncate the covariance beyond 10 times the smoothing distance
            trunc_vertical_smoothing_distance = 10.0 * vertical_smoothing_distance

        if trunc_horizontal_smoothing_distance is None or trunc_horizontal_smoothing_distance <= 0.:
            # truncate the covariance beyond 10 times the smoothing distance
            trunc_horizontal_smoothing_distance = 10. * horizontal_smoothing_distance

        vsd = vertical_smoothing_distance
        hsd = horizontal_smoothing_distance
        tvsd = trunc_vertical_smoothing_distance
        thsd = trunc_horizontal_smoothing_distance

        if vsd == 0.:
            # trick (not optimal) to avoid vertical smoothing I set the truncation distance to 0
            vsd, tvsd = 1.0, 0.

        if hsd == 0.:
            # trick (not optimal) to avoid horizontal smoothing I set the truncation distance to 0
            hsd, thsd = 1.0, 0.

        assert hsd > 0 and vsd > 0
        # ========= load vs uncertainty at each node
        vs_uncertainties = []
        for nnode in range(len(self)):
            vs84_n = depthmodel_from_mod96(self.p84files[nnode]).vs
            vs16_n = depthmodel_from_mod96(self.p16files[nnode]).vs
            vs84_n = depthmodel1D(ztop, vs84_n.interp(zmid))
            vs16_n = depthmodel1D(ztop, vs16_n.interp(zmid))
            vs_unc_n = 0.5 * (vs84_n.values - vs16_n.values)
            vs_unc_n = scale_uncertainties * vs_unc_n + add_uncertainty

            if lock_half_space:
                vs_unc_n[-1] = lock_half_space_sigma

            # make sure all uncertainties are positive
            vs_unc_n = np.clip(vs_unc_n, lock_half_space_sigma, np.inf)
            vs_uncertainties.append(vs_unc_n)

        # ========== compute CM
        CM_row_ind = np.array([], int)
        CM_col_ind = np.array([], int)
        CM_data = np.array([], float)

        for nnode in range(len(self)):
            vs_unc_n = vs_uncertainties[nnode]

            for mnode in range(nnode, len(self)):
                vs_unc_m = vs_uncertainties[mnode]

                # compute the horizontal distance (scalar)
                horizontal_distance_nm = haversine(
                    loni=self.lons[nnode],
                    lati=self.lats[nnode],
                    lonj=self.lons[mnode],
                    latj=self.lats[mnode])

                if horizontal_distance_nm > thsd:
                    # the two cells are farther than trunc_horizontal_smoothing_distance
                    # leave covariance = 0
                    continue

                # compute the vertical distance (2d)
                dz = np.abs(zmid - zmid[:, np.newaxis])
                Itrunc = dz.flat[:] <= tvsd

                if norm == "L2":
                    # square distance
                    d2 = (dz / vsd) ** 2. + (horizontal_distance_nm / hsd) ** 2.

                    # linear correlation coefficient
                    rhonm = np.exp(-0.5 * d2)

                elif norm == "L1":
                    # distance
                    d1 = np.abs(dz / vsd) + np.abs(horizontal_distance_nm / hsd)

                    # linear correlation coefficient
                    rhonm = np.exp(-d1)

                else:
                    raise

                # covariance matrix
                covnm = vs_unc_n[:, np.newaxis] * vs_unc_m * rhonm

                # place the matrix in a large matrix
                col_ind_small, row_ind_small = np.meshgrid(
                    range(len(ztop)), range(len(ztop)))

                row_ind_big = nnode * len(ztop) + row_ind_small.flat[Itrunc]
                col_ind_big = mnode * len(ztop) + col_ind_small.flat[Itrunc]
                covnm = covnm.flat[Itrunc]

                # fill the upper triangle
                CM_row_ind = np.hstack((CM_row_ind, row_ind_big))
                CM_col_ind = np.hstack((CM_col_ind, col_ind_big))
                CM_data = np.hstack((CM_data, covnm))

                if mnode > nnode:
                    # fill the lower triangle
                    CM_row_ind = np.hstack((CM_row_ind, col_ind_big))
                    CM_col_ind = np.hstack((CM_col_ind, row_ind_big))
                    CM_data = np.hstack((CM_data, covnm))

        CM = csc_matrix((CM_data, (CM_row_ind, CM_col_ind)),
                shape=(len(self) * len(ztop), len(self) * len(ztop)))

        # ============== quality control
        # verifies symmetry
        check = (CM - CM.T)
        assert not len(check.data)

        if visual_qc:
            # warning memory error
            # visualize the matrix
            CM_ = CM.toarray()
            CM_ = np.ma.masked_where(CM_ == 0, CM_)
            plt.figure()
            plt.colorbar(plt.imshow(np.log(CM_)))
            plt.show()

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

    def inv_to_mod96strings(self, Model):
        m96strings = []
        for nnode, model in enumerate(self.split(Model)):
            m96strings.append(self.parameterizers[nnode].inv_to_mod96string(model))
        return m96strings

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

    def get_CD(self, scale_uncertainties=1., add_uncertainty=0.):
        CDinv_diag = []
        for n in range(len(self.datacoders)):
            self.datacoders[n].dvalues = scale_uncertainties * self.datacoders[n].dvalues + add_uncertainty

            assert np.all(self.datacoders[n].dvalues > 0.)
            dobs_node, CDinv_diag_node = self.datacoders[n].target()
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
            data = Data[i_current: i_current + ndatapoints]
            i_current += ndatapoints
            datas.append(data)

        assert i_current == len(Data)
        return datas

    def inv(self, Data):
        Values = []
        for nnode, (datacoder, data) in enumerate(zip(self.datacoders, self.split(Data))):
            Values.append(datacoder.inv(data))
        return Values

    def inv_to_laws(self, Data):
        Laws = []
        for nnode, (datacoder, data) in enumerate(zip(self.datacoders, self.split(Data))):
            Laws.append(datacoder.inv_to_laws(data))
        return Laws

    def inv_to_surf96strings(self, Data):
        surf96strings = []
        for nnode, (datacoder, data) in enumerate(zip(self.datacoders, self.split(Data))):
            s96 = datacoder.inv_to_surf96string(data)
            surf96strings.append(s96)
        return surf96strings

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


def sizeformat(size):
    n = 0
    bkmg = []
    while size and n < 4:
        bkmg.append(int(size % 1024))
        size = size // 1024
        n += 1
    return str(bkmg[-1]) + " KMG"[len(bkmg) - 1]


def save_matrix(filename, matrix, verbose):
    if not issparse(matrix):
        sparsity = 0.
        size = matrix.size * matrix.itemsize

    else:
        n = matrix.count_nonzero()
        s = np.prod(matrix.shape)
        sparsity = 100. * (1. - n / float(s))
        size = s * matrix.data.itemsize

    if verbose:
        print('saving {:<30s} {:<10s} {:<14s} {:<10s} {:<10s} {:<10s}'.format(
            filename,
            matrix.__class__.__name__,
            "sparse[{:5.1f}%]".format(sparsity) if sparsity else "dense",
            "x".join(np.asarray(matrix.shape, str)),
            matrix.dtype,
            sizeformat(size)))

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


def chi2_data(niter=None, Data=None, Dobs=None, CDinv=None):
    if Data is None:
        Data = np.load(DFILE.format(niter=niter))
    if Dobs is None:
        Dobs = np.load(DOBSFILE)
    if CDinv is None:
        CDinv = load_sparse_npz(CDINVFILE.format(niter))
    assert issparse(CDinv)
    assert not issparse(Data - Dobs)

    Inan = np.isnan(Data) | np.isnan(Dobs)
    Dobs[Inan] = Data[Inan] = 0.  # so that they do not act on the data misfit

    return np.dot((Data - Dobs), CDinv * (Data - Dobs))  # only because CDinv is known exactly


def chi2_model(niter=None, Model=None, Mprior=None, CM=None):

    assert issparse(CM)
    assert not issparse(Model - Mprior)

    if Model is None:
        Model = np.load(MFILE.format(niter=niter))
    if Mprior is None:
        Mprior = np.load(MPRIORFILE)
    if CM is None:
        CM = load_sparse_npz(CMFILE.format(niter))

    return np.dot(Model - Mprior, Ainv_dot_b(A=CM, b=Model - Mprior))


def print_costs(niter, data_cost=None, model_cost=None):
    if data_cost is None:
        data_cost = chi2_data(niter=niter)
    if model_cost is None:
        model_cost = chi2_model(niter=niter)

    print('iter {:03d}: chi2_data={:<16f} chi2_model={:<16f} chi2={:<16f}'.format(
        niter, data_cost, model_cost, data_cost + model_cost))


def load_last_frechet_derivatives(supertheory, verbose, mapkwargs):
    niter = lastiter()
    while niter >= 0:
        fdfilename = FDFILE.format(niter=niter)
        if os.path.isfile(fdfilename):
            if verbose:
                print('loading {}'.format(fdfilename))
            Gi = load_sparse_npz(fdfilename)
            return Gi
        niter -= 1

    if niter < 0:
        # no frechet derivatives found even G0,
        return update_frechet_derivatives(supertheory, verbose, mapkwargs)


def update_frechet_derivatives(supertheory, verbose, mapkwargs):
    """

    :param supertheory:
    :param verbose:
    :param mapkwargs:
    :return:
    """
    # if os.path.isfile(FDFILE.format(niter=0)):
    #     # test do not update the FD
    #     G0 = load_sparse_npz(FDFILE.format(niter=0))
    #     return G0

    # ========= compute the frechet derivatives for the last model found
    niter = lastiter()
    mfilename = MFILE.format(niter=niter)
    fdfilename = FDFILE.format(niter=niter)
    if not os.path.isfile(fdfilename):
        if verbose:
            print('computing frechet derivatives for iteration {}...'.format(niter))
        Mi = np.load(mfilename)
        Gi = supertheory.get_FD(Mi, mapkwargs)
        save_matrix(fdfilename, Gi, verbose)
        return Gi

    else:
        if verbose:
            print('{} exists already'.format(fdfilename))
        Gi = load_sparse_npz(fdfilename)
        return Gi


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
            m96fileout = M96FILEOUT.format(node="*", niter=niter)
            s96fileout = S96FILEOUT.format(node="*", niter=niter)
            trashcmd = 'trash {} {} {} {} {}'.format(mfile, dfile, fdfile, m96fileout, s96fileout)
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
        horizontal_smoothing_distance = argv["-prior"][0]
        vertical_smoothing_distance = argv["-prior"][1]
        trunc_horizontal_smoothing_distance = argv["-prior"][2]
        trunc_vertical_smoothing_distance = argv["-prior"][3]
        lock_half_space = bool(int(argv["-prior"][4]))
        scale_uncertainties = argv["-prior"][5]
        add_uncertainty = argv["-prior"][6]

        # ========= set the prior model and covariance matrix
        Mprior = superparameterizer.get_Mprior()
        CM = nodefile.get_CM(
            horizontal_smoothing_distance=horizontal_smoothing_distance,
            vertical_smoothing_distance=vertical_smoothing_distance,
            trunc_horizontal_smoothing_distance=trunc_horizontal_smoothing_distance,
            trunc_vertical_smoothing_distance=trunc_vertical_smoothing_distance,
            lock_half_space=lock_half_space,
            scale_uncertainties=scale_uncertainties,
            add_uncertainty=add_uncertainty,
            norm="L1",
            visual_qc=False)

        save_matrix(MPRIORFILE, Mprior, verbose)
        save_matrix(CMFILE, CM, verbose)

        M0 = Mprior
        m0filename = MFILE.format(niter=0)
        save_matrix(m0filename, M0, verbose)

        d0filename = DFILE.format(niter=0)
        D0 = supertheory(M0)
        save_matrix(d0filename, D0, verbose)

    if "-data" in argv.keys():
        if len(argv['-data']) == 2:
            scale_uncertainties, add_uncertainty = argv["-data"]
        elif len(argv['-data']) == 0:
            scale_uncertainties, add_uncertainty = 1.0, 0.
        else:
            raise ValueError('see help for -data')

        # ========= set the observed data and data covariance matrix
        Dobs = superdatacoder.get_Dobs()

        CD, CDinv = superdatacoder.get_CD(
            scale_uncertainties=scale_uncertainties,
            add_uncertainty=add_uncertainty)

        save_matrix(DOBSFILE, Dobs, verbose)
        save_matrix(CDFILE, CD, verbose)
        save_matrix(CDINVFILE, CDinv, verbose)

    if "-fd" in argv.keys():
        # manually update the frechet derivatives
        update_frechet_derivatives(supertheory, verbose, mapkwargs)

    if "-upd" in argv.keys():
        number_of_iterations = argv["-upd"][0] if len(argv["-upd"]) > 0 else 1
        update_G = bool(argv["-upd"][1]) if len(argv["-upd"]) > 1 else True

        # ================ load data
        Mprior = np.load(MPRIORFILE)
        M0 = np.load(MFILE.format(niter=0))
        Dobs = np.load(DOBSFILE)

        CM = load_sparse_npz(CMFILE)
        CD = load_sparse_npz(CDFILE)
        CDinv = load_sparse_npz(CDINVFILE)

        for _ in range(number_of_iterations):
            niter = lastiter()

            if update_G:
                # compute frechet derivatives for this iteration (if not already computed)
                Gi = update_frechet_derivatives(supertheory, verbose=verbose, mapkwargs=mapkwargs)
            else:
                # find the last version of the frechet derivatives, compute G0 if none is found
                Gi = load_last_frechet_derivatives(supertheory, verbose=verbose, mapkwargs=mapkwargs)

            # ==== save the current model state
            mfilename = MFILE.format(niter=niter)
            dfile = DFILE.format(niter=niter)
            # fdfilename = FDFILE.format(niter=niter)

            Mi = np.load(mfilename)
            Di = np.load(dfile)  # g(Mi)
            # Gi = load_sparse_npz(fdfilename)

            # ==== update the current model
            has_nan = np.isnan(Dobs).any()
            Dinew, Minew = tv23_1(Dobs, Di, CD,
                       Mprior, M0, Mi, CM,
                       Gi, supertheory)
            if has_nan:
                assert np.isnan(Dobs).any()

            # ==== save the new model state
            mfilename_new = MFILE.format(niter=niter + 1)
            dfilename_new = DFILE.format(niter=niter + 1)
            # fdfilename_new = FDFILE.format(niter=niter + 1)

            save_matrix(mfilename_new, Minew, verbose)
            save_matrix(dfilename_new, Dinew, verbose)
            # save_matrix(fdfilename_new, Ginew, verbose)

            if verbose:
                data_cost = chi2_data(niter=None, Data=Dinew, Dobs=Dobs, CDinv=CDinv)
                model_cost = chi2_model(niter=None, Model=Minew, Mprior=Mprior, CM=CM)
                print_costs(niter+1, data_cost=data_cost, model_cost=model_cost)

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

            data_cost = chi2_data(None, Data=Data, Dobs=Dobs, CDinv=CDinv)
            model_cost = chi2_model(None, Model=Model, Mprior=Mprior, CM=CM)

            if verbose:
                print_costs(niter, data_cost=data_cost, model_cost=model_cost)

            # model_cost = np.dot(np.dot((Model - Mprior), CMinv), (Model - Mprior))  # NO
            data_costs.append(data_cost)
            model_costs.append(model_cost)

        ax_costs1.plot(data_costs, 'ko')
        ax_costs2.plot(model_costs, 'bo')
        ax_costs3.plot(np.array(data_costs) + np.array(model_costs), 'bo')

        # =================== display the models side by side
        fig_model = plt.figure()
        ax_model = fig_model.add_subplot(111)

        try:
            superparameterizer.show_models(ax_model, Mprior - Mpriorunc, color="k", alpha=0.3)
            superparameterizer.show_models(ax_model, Mprior + Mpriorunc, color="k", alpha=0.3)
            superparameterizer.show_models(ax_model, Mprior, Munc=Mpriorunc, color="k", alpha=0.3)
        except ValueError as e:
            warnings.warn(str(e))
        superparameterizer.show_models(ax_model, Mprior, color="k", alpha=1.0)

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

    if "-save" in argv.keys():
        niter = lastiter()
        Model = np.load(MFILE.format(niter=niter))
        Data = np.load(DFILE.format(niter=niter))
        m96strings = superparameterizer.inv_to_mod96strings(Model)
        s96strings = superdatacoder.inv_to_surf96strings(Data)
        for nnode, (m96string, s96string) in enumerate(zip(m96strings, s96strings)):
            node_name = nodefile.nodes[nnode]

            m96fileout = M96FILEOUT.format(node=node_name, niter=niter)
            s96fileout = S96FILEOUT.format(node=node_name, niter=niter)

            if verbose:
                print('writing {}'.format(m96fileout))
                print('writing {}'.format(s96fileout))

            with open(m96fileout, 'w') as fid:
                fid.write(m96string)

            with open(s96fileout, 'w') as fid:
                fid.write(s96string)


def tv23_1(Dobs, Di, CD,
           Mprior, M0, Mi, CM,
           Gi, supertheory):

    # nan issues
    # nan can occur in the predicted data if a data point is above the cut off period
    # let ignore them
    Inan = np.isnan(Di) | np.isnan(Dobs)
    Dobs[Inan] = Di[Inan] = 0.
    Xi = Dobs - Di + Gi * (Mi - Mprior)

    CMGiT = CM * Gi.T

    Ai = CD + Gi * CMGiT

    if False:
        Aiinv_dot_Xi = spsolve(A=Ai, b=Xi)
    else:
        lu = splu(Ai)
        Aiinv_dot_Xi = lu.solve(Xi)

    error = np.abs(Xi - Ai * Aiinv_dot_Xi).sum()
    print('error on A.x=b : {}'.format(error))

    KiXi = CMGiT * Aiinv_dot_Xi

    mu = 1.0
    Minew = Dinew = None
    while mu > 0.001:
        try:
            Minew = Mi + mu * (M0 + KiXi - Mi)
            Dinew = supertheory(Minew)
            break
        except CPiSDomainError as e:
            print('the step was too large and leaded to an error in the forward problem \n({}), '
                  'try reducing the step length'.format(str(e)))
            mu /= 2.0
            continue

    if Minew is None:
        raise Exception('fail')

    return Dinew, Minew
