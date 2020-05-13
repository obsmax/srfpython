from __future__ import print_function
from builtins import input
import warnings
import glob, os
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import spsolve, splu  #inv as sparse_inv NO!!!
import matplotlib.pyplot as plt
from srfpython.HerrMet.nodefile import NodeFile
from srfpython.depthdisp.depthmodels import depthmodel_from_mod96, depthmodel1D
from srfpython.HerrMet.paramfile import load_paramfile
from srfpython.HerrMet.datacoders import makedatacoder, Datacoder_log
from srfpython.HerrMet.theory import Theory
from srfpython.standalone.multipro8 import Job, MapSync
from srfpython.standalone.stdout import waitbarpipe
from srfpython.Herrmann.Herrmann import CPiSDomainError
from srfpython.HerrMet.files import ROOTKEY
from srfpython.coordinates import haversine
import pickle

# TODO add convergence interruption warning
# TODO speed up CM computations
# TODO add safety gards about the size of CM

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
CMFILE = os.path.join(WORKINGDIR, ROOTKEY + "CM_triu.npz")

SUPERPARAMETERIZERFILE = os.path.join(WORKINGDIR, ROOTKEY + "SP.pkl")
SUPERDATACODERFILE = os.path.join(WORKINGDIR, ROOTKEY + "SD.pkl")
SUPERTHOERYFILE = os.path.join(WORKINGDIR, ROOTKEY + "ST.pkl")

M96FILEOUT = os.path.join(WORKINGDIR, ROOTKEY + "{node}_optimized.mod96")
S96FILEOUT = os.path.join(WORKINGDIR, ROOTKEY + "{node}_optimized.surf96")
M96FILEOUTS = os.path.join(WORKINGDIR, ROOTKEY + "*_optimized.mod96")
S96FILEOUTS = os.path.join(WORKINGDIR, ROOTKEY + "*_optimized.surf96")

CLEAR_LIST = [
     NODEFILELOCAL,
     MPRIORFILE,
     DOBSFILE,
     CDFILE,
     CDINVFILE,
     CMFILE,
     MFILES,
     DFILES,
     FDFILES,
     M96FILEOUTS,
     S96FILEOUTS,
     SUPERPARAMETERIZERFILE,
     SUPERDATACODERFILE,
     SUPERTHOERYFILE]


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
authorized_keys = ["-init", "-restart", "-prior",
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
    -prior 80. 1. 160. 2. 1 1.0 0.0   \\
    -data 1.0 0.0
    
# the prior can be changed afterwards, 
# it requires to remove the iteration files using option -restart
HerrMet --optimize \\
    -prior 100. 2. 200. 4. 1 1.0 0.0    

# run 1 iterations with fd update
HerrMet --optimize -upd 1 1

# run 2 iterations without fd update
HerrMet --optimize -upd 2 0

# update the fd manually
HerrMet --optimize -fd 
    
# run 1 more iteration without fd update 
HerrMet --optimize -upd 1 0   # -upd 1 1 would be equivalent

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

    def get_superparameterizer(self, verbose, mapkwargs):
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

        def job_generator():
            for nnode, (node, lon, lat, targetfile, paramfile, medianfile, p16file, p84file) in enumerate(self):
                yield Job(medianfile, self.ztop, self.zmid)

        def job_handler(medianfile, ztop, zmid):
            dm = depthmodel_from_mod96(medianfile)
            vs = dm.vs.interp(zmid)

            parameter_string = parameter_string_header

            for i in range(len(ztop)):
                # SET VINF < VS extracted from pointwise inv < VSUP
                # such as parameterizer.MMEAN corresponds to the extracted vs
                parameter_string += "VS{} {} {}\n".format(i, vs[i]-0.01, vs[i]+0.01)
            parameterizer = load_paramfile(parameter_string, verbose=False)[0]
            return parameter_string, parameterizer

        wb = None
        if verbose:
            wb = waitbarpipe('parameterizers')

        parameterizers = []
        with MapSync(job_handler, job_generator(), **mapkwargs) as ma:
            for jobid, (parameter_string, parameterizer), _, _ in ma:
                parameterizers.append(parameterizer)

                if verbose:
                    wb.refresh(jobid / float(len(self)))

        if verbose:
            wb.close()

        return SuperParameterizer(parameterizers)

    def get_superdatacoder(self, verbose):
        """

        :param self:
        :return:
        """
        datacoder_strings = []

        if verbose:
            wb = waitbarpipe('datacoders    ')

        for nnode, (node, lon, lat, targetfile, paramfile, medianfile, p16file, p84file) in enumerate(self):

            with open(targetfile, 'r') as fid:
                datacoder_string = "".join(fid.readlines())
                datacoder_strings.append(datacoder_string)

            if verbose:
                wb.refresh(nnode / float(len(self)))

        if verbose:
            wb.close()

        datacoder_strings = np.asarray(datacoder_strings, str)
        datacoders = [makedatacoder(datacoder_string, which=Datacoder_log) for datacoder_string in datacoder_strings]
        return SuperDatacoder(datacoders)

    def build_or_load_forward_problem(self, verbose, mapkwargs):
        """

        :param verbose:
        :param mapkwargs:
        :return:
        """
        if verbose:
            print('building the forward problem operator ...')

        if os.path.isfile(SUPERPARAMETERIZERFILE):
            if verbose:
                print('loading ', SUPERPARAMETERIZERFILE)
            with open(SUPERPARAMETERIZERFILE, 'rb') as fid:
                superparameterizer = pickle.load(fid)

        else:
            superparameterizer = self.get_superparameterizer(verbose, mapkwargs)
            with open(SUPERPARAMETERIZERFILE, 'wb') as fid:
                print('saving ', SUPERPARAMETERIZERFILE)
                fid.write(pickle.dumps(superparameterizer))

        if os.path.isfile(SUPERDATACODERFILE):
            if verbose:
                print('loading ', SUPERDATACODERFILE)
            with open(SUPERDATACODERFILE, 'rb') as fid:
                superdatacoder = pickle.load(fid)

        else:
            superdatacoder = self.get_superdatacoder(verbose)
            with open(SUPERDATACODERFILE, 'wb') as fid:
                print('saving ', SUPERDATACODERFILE)
                fid.write(pickle.dumps(superdatacoder))

        if os.path.isfile(SUPERTHOERYFILE):
            if verbose:
                print('loading ', SUPERTHOERYFILE)
            with open(SUPERTHOERYFILE, 'rb') as fid:
                supertheory = pickle.load(fid)
            supertheory.mapkwargs = mapkwargs  # set the new mapkwargs to supertheory

        else:
            supertheory = SuperTheory(
                superparameterizer=superparameterizer,
                superdatacoder=superdatacoder,
                verbose=verbose,
                mapkwargs=mapkwargs)
            with open(SUPERTHOERYFILE, 'wb') as fid:
                print('saving ', SUPERTHOERYFILE)
                fid.write(pickle.dumps(supertheory))

        n_data_points, n_parameters = supertheory.shape

        if verbose:
            print('done')

        return n_data_points, n_parameters, superdatacoder, superparameterizer, supertheory

    def get_CM(
            self, horizontal_smoothing_distance, vertical_smoothing_distance,
            trunc_horizontal_smoothing_distance=0.,
            trunc_vertical_smoothing_distance=0.,
            lock_half_space=True,
            lock_half_space_sigma=0.01,
            scale_uncertainties=1.,
            add_uncertainty=0.,
            norm="L1",
            visual_qc=False,
            verbose=True,
            mapkwargs=None):

        if norm != "L1":
            raise NotImplementedError('L2 not recommended (unstable)')
        vsd = vertical_smoothing_distance
        hsd = horizontal_smoothing_distance
        tvsd = trunc_vertical_smoothing_distance
        thsd = trunc_horizontal_smoothing_distance

        if not vsd >= 0 and tvsd >= 0 and hsd >= 0 and thsd >= 0:
            raise ValueError(vsd, tvsd, hsd, thsd)

        # ========= load vs uncertainty array at each node
        print('    get_vs_uncertainties')
        vs_uncertainties = self.get_vs_uncertainties(
            scale_uncertainties=scale_uncertainties,
            add_uncertainty=add_uncertainty,
            lock_half_space=lock_half_space,
            lock_half_space_sigma=lock_half_space_sigma,
            mapkwargs=mapkwargs)
        print('    done')

        # =========== Compute CM_triu with truncature
        nnodes = len(self)
        nlayer = len(self.ztop)
        n_parameters = nlayer * nnodes
        vsmooth = vsd > 0 and tvsd > 0
        hsmooth = hsd > 0 and thsd > 0
        if not (vsmooth or hsmooth):
            # no smoothing at all, easy
            CM_triu = sp.diags(vs_uncertainties.flat[:] ** 2.0,
                       shape=(n_parameters, n_parameters),
                       format="csc", dtype=float)
        else:

            # =========== prepare smoothing
            print('    prepare_vertical_smoothing')
            rhoz_nupper, rhoz_nlower, rhoz_triu, rhoz_tril = \
                self.prepare_vertical_smoothing(vsd, tvsd, norm=norm)
            print('    done')

            print('    prepare_horizontal_smoothing')
            rhod_triu = self.prepare_horizontal_smoothing(hsd, thsd, norm=norm, mapkwargs=mapkwargs)
            print('    done')

            # ===========
            CM_triu_rows = []
            CM_triu_cols = []
            CM_triu_data = []

            if verbose:
                wb = waitbarpipe()

            for nnode in range(len(self)):
                vs_unc_n = vs_uncertainties[nnode, :]

                if vsmooth:

                    covnn_triu_row = nnode * nlayer + rhoz_triu.row
                    covnn_triu_col = nnode * nlayer + rhoz_triu.col
                    covnn_triu_data = vs_unc_n[rhoz_triu.row] * vs_unc_n[rhoz_triu.col] * rhoz_triu.data

                    CM_triu_rows.append(covnn_triu_row)
                    CM_triu_cols.append(covnn_triu_col)
                    CM_triu_data.append(covnn_triu_data)

                    if hsmooth:

                        for mnode in range(nnode + 1, len(self)):
                            if rhod_triu[nnode, mnode] == 0.:
                                continue

                            vs_unc_m = vs_uncertainties[mnode, :]

                            if 1:
                                covnm_row = np.zeros(rhoz_nupper + rhoz_nlower, int)
                                covnm_col = np.zeros(rhoz_nupper + rhoz_nlower, int)
                                covnm_data = np.zeros(rhoz_nupper + rhoz_nlower, float)

                                covnm_row[:rhoz_nupper] = nnode * nlayer + rhoz_triu.row
                                covnm_col[:rhoz_nupper] = mnode * nlayer + rhoz_triu.col
                                covnm_data[:rhoz_nupper] = vs_unc_n[rhoz_triu.row] * vs_unc_m[rhoz_triu.col] \
                                                  * rhoz_triu.data * rhod_triu[nnode, mnode]

                                covnm_row[rhoz_nupper:] = nnode * nlayer + rhoz_tril.row
                                covnm_col[rhoz_nupper:] = mnode * nlayer + rhoz_tril.col
                                covnm_data[rhoz_nupper:] = vs_unc_n[rhoz_tril.row] * vs_unc_m[rhoz_tril.col] \
                                                  * rhoz_tril.data * rhod_triu[nnode, mnode]

                                CM_triu_rows.append(covnm_row)
                                CM_triu_cols.append(covnm_col)
                                CM_triu_data.append(covnm_data)
                            else:
                                # only accounts for the diagonal terms since nnode != mnode
                                raise Exception('seems wrong : inversion fails')
                                covnm_row = nnode * nlayer + np.arange(nlayer)
                                covnm_col = mnode * nlayer + np.arange(nlayer)
                                covnm_data = vs_unc_n * vs_unc_m * rhod_triu[nnode, mnode]
                                CM_triu_rows.append(covnm_row)
                                CM_triu_cols.append(covnm_col)
                                CM_triu_data.append(covnm_data)

                else:
                    assert hsmooth  # implicit
                    covnn_row = nnode * nlayer + np.arange(nlayer)
                    covnn_col = nnode * nlayer + np.arange(nlayer)
                    covnn_data = vs_unc_n ** 2.

                    CM_triu_rows.append(covnn_row)
                    CM_triu_cols.append(covnn_col)
                    CM_triu_data.append(covnn_data)

                    for mnode in range(nnode + 1, len(self)):
                        vs_unc_m = vs_uncertainties[mnode, :]

                        covnm_row = nnode * nlayer + np.arange(nlayer)
                        covnm_col = mnode * nlayer + np.arange(nlayer)
                        covnm_data = vs_unc_n * vs_unc_m * rhod_triu[nnode, mnode]

                        CM_triu_rows.append(covnm_row)
                        CM_triu_cols.append(covnm_col)
                        CM_triu_data.append(covnm_data)

                if verbose:
                    wb.refresh(nnode / float(len(self)))

            if verbose:
                wb.close()
            CM_triu_rows = np.hstack(CM_triu_rows)
            CM_triu_cols = np.hstack(CM_triu_cols)
            CM_triu_data = np.hstack(CM_triu_data)
            CM_triu = sp.csc_matrix((CM_triu_data, (CM_triu_rows, CM_triu_cols)),
                           shape=(n_parameters, n_parameters), dtype=float)

        # =========== fix truncature issues,
        # approximate CM = Vp * Lp * Vp.T
        # where Vp are the first eigenvectors
        #       Lp is strictly diagonal and positive
        #       Vp.T * Vp = Id
        #       Vp * Vp.T != Id
        from scipy.sparse.linalg import eigsh

        CM = CM_triu + sp.tril(CM_triu.T, k=-1, format="csc")  # add lower triangle without diag (already in triu)
        CM_first_eigenvalues, CM_first_eigenvectors = eigsh(CM, k=CM.shape[0] // 4)
        del CM
        I = CM_first_eigenvalues > 0.

        CM_Vp = sp.csc_matrix(CM_first_eigenvectors[:, I])
        CM_Lp = sp.diags(CM_first_eigenvalues[I], format="csc")

        CM_triu = sp.triu(CM_Vp * CM_Lp * CM_Vp.T)

        if visual_qc:
            CM = CM_triu + sp.tril(CM_triu.T, k=-1, format="csc")  # add lower triangle without diag (already in triu)

            if CM.shape[0] > 2000:
                raise Exception('CM is too large')
            A = CM.A
            plt.figure()
            A = np.ma.masked_where(A == 0, A)
            plt.imshow(A)
            plt.show()

        return CM_triu

    def prepare_horizontal_smoothing(self, hsd, thsd, norm, mapkwargs):
        if hsd == 0 and thsd == 0:
            rhod_triu = None

        else:
            Nnodes = len(self)

            # warning : I do not fill the lower triangle and diagonal !!!
            rhod_rows = [np.arange(Nnodes)]
            rhod_cols = [np.arange(Nnodes)]
            rhod_datas = [np.ones(Nnodes)]

            # fill non-diagonal terms
            def job_generator():
                for nnode in range(Nnodes - 1):
                    yield Job(nnode,
                              loni=self.lons[nnode],
                              lati=self.lats[nnode],
                              lonj=self.lons[nnode + 1:],
                              latj=self.lats[nnode + 1:])

            def job_handler(nnode, **kwargs):
                horizontal_distances = haversine(**kwargs)
                mnodes = np.arange(nnode+1, Nnodes)

                I = horizontal_distances < thsd
                nnodes = nnode * np.ones(I.sum())
                mnodes = mnodes[I]
                horizontal_distances = horizontal_distances[I]
                if norm == "L1":
                    datas = np.exp(-horizontal_distances / hsd)
                elif norm == "L2":
                    datas = np.exp(-horizontal_distances / hsd)
                else:
                    raise NotImplementedError

                return nnodes, mnodes, datas

            with MapSync(job_handler, job_generator(), **mapkwargs) as ma:
                for jobid, (nnodes, mnodes, datas), _, _ in ma:
                    rhod_rows.append(nnodes)
                    rhod_cols.append(mnodes)
                    rhod_datas.append(datas)

            rhod_rows = np.hstack(rhod_rows)
            rhod_cols = np.hstack(rhod_cols)
            rhod_datas = np.hstack(rhod_datas)

            rhod_triu = sp.csc_matrix((rhod_datas, (rhod_rows, rhod_cols)),
                                      shape=(Nnodes, Nnodes), dtype=float)
        return rhod_triu

    def prepare_vertical_smoothing(self, vsd, tvsd, norm):
        if vsd == 0 and tvsd == 0:
            rhoz_nupper = rhoz_nlower = rhoz_triu = rhoz_tril = None
        else:
            zmid = self.zmid
            dz = np.abs(zmid - zmid[:, np.newaxis])
            if norm == "L1":
                rhoz = np.exp(-dz / vsd)
            elif norm == "L2":
                rhoz = np.exp(-(dz / vsd) ** 2.0)
            else:
                raise NotImplementedError
            rhoz[dz > tvsd] = 0.
            rhoz_triu = sp.triu(rhoz, format="coo", k=0)  # upper triangle with diagonal
            rhoz_tril = sp.tril(rhoz, format="coo", k=-1)  # lower triangle without diagonal
            rhoz_nupper, rhoz_nlower = len(rhoz_triu.data), len(rhoz_tril.data)

        return rhoz_nupper, rhoz_nlower, rhoz_triu, rhoz_tril

    def get_vs_uncertainties(self, scale_uncertainties, add_uncertainty, lock_half_space, lock_half_space_sigma, mapkwargs):
        ztop = self.ztop
        zmid = self.zmid
        nnodes = len(self)
        nlayer = len(ztop)

        vs_uncertainties = np.zeros((nnodes, nlayer), float)  # one row per node, one col per layer
        #
        # def job_generator():
        #     for nnode in range(nnodes):
        #         p16file, p84file = self.p16files[nnode], self.p84files[nnode]
        #         yield Job(nnode, p16file, p84file)
        #
        # def job_handler(nnode, p16file, p84file):
        #     vs84_n = depthmodel_from_mod96(p84file).vs
        #     vs16_n = depthmodel_from_mod96(p16file).vs
        #     vs84_n = depthmodel1D(ztop, vs84_n.interp(zmid))
        #     vs16_n = depthmodel1D(ztop, vs16_n.interp(zmid))
        #     vs_unc_n = 0.5 * (vs84_n.values - vs16_n.values)
        #     vs_unc_n = scale_uncertainties * vs_unc_n + add_uncertainty
        #
        #     if lock_half_space:
        #         vs_unc_n[-1] = lock_half_space_sigma
        #
        #     # make sure all uncertainties are positive
        #     vs_unc_n = np.clip(vs_unc_n, lock_half_space_sigma, np.inf)
        #     return nnode, vs_unc_n
        #
        # if mapkwargs is None:
        #     mapkwargs = {}
        # with MapSync(job_handler, job_generator(), **mapkwargs) as ma:
        #     for jobid, (nnode, vs_unc_n), _, _ in ma:
        #         vs_uncertainties[nnode, :] = vs_unc_n
        # return vs_uncertainties


        for nnode in range(nnodes):
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
            vs_uncertainties[nnode, :] = vs_unc_n
        return vs_uncertainties


class SuperParameterizer(object):

    def __init__(self, parameterizers):
        self.parameterizers = parameterizers

    def __len__(self):
        return len(self.parameterizers)

    def get_Mprior(self):
        Mprior = []
        for parameterizer in self.parameterizers:
            mprior = parameterizer.MMEAN
            Mprior.append(mprior)
        return np.hstack(Mprior)

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

    def __len__(self):
        return len(self.datacoders)

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

        CDinv = sp.diags(CDinv_diag, offsets=0, format="csc", dtype=float)
        CD = sp.diags(CDinv_diag ** -1.0, offsets=0, format="csc", dtype=float)
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
    def __init__(self, superparameterizer, superdatacoder, verbose, mapkwargs):
        assert isinstance(superdatacoder, SuperDatacoder)
        assert isinstance(superparameterizer, SuperParameterizer)
        assert len(superparameterizer) == len(superdatacoder)
        self.mapkwargs = mapkwargs

        self.superparameterizer = superparameterizer
        self.superdatacoder = superdatacoder

        def job_generator():
            for nnode, (parameterizer, datacoder) in enumerate(zip(
                    superparameterizer.parameterizers, superdatacoder.datacoders)):
                yield Job(nnode, parameterizer, datacoder)

        def job_handler(nnode, parameterizer, datacoder):
            theory = Theory(parameterizer=parameterizer, datacoder=datacoder)
            return nnode, theory

        if verbose:
            wb = waitbarpipe('forward operators')

        self.theorys = []
        with MapSync(job_handler, job_generator(), **mapkwargs) as ma:
            for jobid, (nnode, theory), _, _ in ma:
                self.theorys.append(theory)

                if verbose:
                    wb.refresh(jobid / float(len(self.superdatacoder)))

        if verbose:
            wb.close()

        # number of data points for each node
        self.nds = np.asarray([len(theory.datacoder.values) for theory in self.theorys], int)

        # number of parameters for each node
        self.nps = np.asarray([len(theory.parameterizer.MMEAN) for theory in self.theorys], int)

    def __len__(self):
        return len(self.superdatacoder)

    @property
    def shape(self):
        return self.nds.sum(), self.nps.sum()

    def __call__(self, Model):
        models = self.superparameterizer.split(Model)

        def job_generator():
            for theory, model in zip(self.theorys, models):
                yield Job(theory, model)

        def job_handler(theory, model):
            data = theory(model)
            return data

        Data = []
        with MapSync(job_handler, job_generator(), **self.mapkwargs) as ma:
            for jobid, data, _, _ in ma:
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

        G = sp.csc_matrix((G_fd_data, (G_row_ind, G_col_ind)), shape=G_shape)
        return G

    def get_CMGiT(self, CM, Gi, mapkwargs, verbose):
        if verbose:
            print('computing CM . G.T ...')
            
        assert Gi.shape == self.shape

        d_end = np.cumsum(self.nds)  # excluded
        d_begin = np.hstack((0, d_end[:-1]))  # included
        p_end = np.cumsum(self.nps)  # excluded
        p_begin = np.hstack((0, p_end[:-1]))  # included
        assert len(p_begin) == len(p_end) == len(d_begin) == len(d_end)
        assert Gi.shape == (d_end[-1], p_end[-1])

        def job_generator():
            for mnode, (pb, pe, db, de) in enumerate(zip(p_begin, p_end, d_begin, d_end)):
                GimmT = Gi[db:de, pb:pe].T
                CM_supercol = CM[:, pb:pe]
                yield Job(pb, pe, db, de, GimmT, CM_supercol)

        def job_handler(pb, pe, db, de, GimmT, CM_supercol):

            rows = []
            cols = []
            datas = []

            for qb, qe in zip(p_begin, p_end):
                # CMnm = CM[qb:qe, pb:pe]
                CMGiTnm = (CM_supercol[qb:qe, :] * GimmT).tocoo()
                if len(CMGiTnm.row):
                    rows.append(qb + CMGiTnm.row)
                    cols.append(db + CMGiTnm.col)
                    datas.append(CMGiTnm.data)

            rows = np.hstack(rows)
            cols = np.hstack(cols)
            datas = np.hstack(datas)

            return rows, cols, datas

        CMGiT_rows = []
        CMGiT_cols = []
        CMGiT_data = []
        if verbose:
            wb = waitbarpipe('')

        with MapSync(job_handler, job_generator(), **mapkwargs) as ma:
            for jobid, (rows, cols, datas), _, _ in ma:

                if len(rows):
                    CMGiT_rows.append(rows)
                    CMGiT_cols.append(cols)
                    CMGiT_data.append(datas)
                if verbose:
                    wb.refresh(jobid / float(len(p_begin)))
        if verbose:
            wb.close()

        CMGiT_rows = np.hstack(CMGiT_rows)
        CMGiT_cols = np.hstack(CMGiT_cols)
        CMGiT_data = np.hstack(CMGiT_data)

        CMGiT = sp.csc_matrix((CMGiT_data, (CMGiT_rows, CMGiT_cols)),
                              shape=(CM.shape[0], Gi.shape[0]),
                              dtype=float)
        if False:
            assert input('sure?') == "y"
            assert np.all(CMGiT.toarray() == (CM * Gi.T).toarray())
            plt.figure()
            plt.subplot(121)
            plt.imshow(CMGiT.toarray())
            plt.subplot(122)
            plt.imshow((CM * Gi.T).toarray())
            plt.show()
        print('ok')
        return CMGiT
        #
        #
        #
        # CMGiT_rows = []
        # CMGiT_cols = []
        # CMGiT_data = []
        #
        # for mnode, (pb, pe, db, de) in enumerate(zip(p_begin, p_end, d_begin, d_end)):
        #     # the "super colonne" is CM[:, jb:je]
        #     GimmT = Gi[db:de, pb:pe].T
        #
        #     for nnode, (qb, qe) in enumerate(zip(p_begin, p_end)):
        #         CMnm = CM[qb:qe, pb:pe]
        #         if CMnm.count_nonzero() == 0:
        #             continue
        #
        #         CMGiTnm = (CMnm * GimmT).tocoo()
        #
        #         if len(CMGiTnm.row):
        #             CMGiT_rows.append(qb + CMGiTnm.row)
        #             CMGiT_cols.append(db + CMGiTnm.col)
        #             CMGiT_data.append(CMGiTnm.data)
        #
        # CMGiT_rows = np.hstack(CMGiT_rows)
        # CMGiT_cols = np.hstack(CMGiT_cols)
        # CMGiT_data = np.hstack(CMGiT_data)
        #
        # CMGiT = sp.csc_matrix((CMGiT_data, (CMGiT_rows, CMGiT_cols)),
        #                       shape=(CM.shape[0], Gi.shape[0]),
        #                       dtype=float)
        # if False:
        #     assert input('sure?') == "y"
        #     assert np.all(CMGiT.toarray() == (CM * Gi.T).toarray())
        #     plt.figure()
        #     plt.subplot(121)
        #     plt.imshow(CMGiT.toarray())
        #     plt.subplot(122)
        #     plt.imshow((CM * Gi.T).toarray())
        #     plt.show()
        # print('ok')
        # return CMGiT


def load_CM(CMFILE):
    # load CM_triu and compute full CM 

    CM_triu = sp.load_npz(CMFILE)
    if sp.isspmatrix_dia(CM_triu):
        CM = CM_triu

    else:
        CM = CM_triu + CM_triu.T - sp.diags(CM_triu.diagonal())

    return CM

#
# def sparsedet(M):
#     lu = splu(M)
#     diagL = lu.L.diagonal()
#     diagU = lu.U.diagonal()
#     diagL = diagL.astype(np.complex128)
#     diagU = diagU.astype(np.complex128)
#     logdet = np.log(diagL).sum() + np.log(diagU).sum()
#     det = np.exp(logdet)  # usually underflows/overflows for large matrices
#     return det.real


def sizeformat(size):
    n = 0
    bkmg = []
    while size and n < 4:
        bkmg.append(int(size % 1024))
        size = size // 1024
        n += 1
    return str(bkmg[-1]) + " KMG"[len(bkmg) - 1]


def save_matrix(filename, matrix, verbose):
    if not sp.issparse(matrix):
        sparsity = 0.
        size = matrix.size * matrix.itemsize

    else:
        n = matrix.count_nonzero()
        s = np.prod(matrix.shape)
        sparsity = 100. * (1. - n / float(s))
        size = n * matrix.data.itemsize

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

    elif isinstance(matrix, sp.spmatrix):
        if not filename.endswith('.npz'):
            raise ValueError('{} does not end with .npz'.format(filename))
        sp.save_npz(filename, matrix)

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
    if sp.issparse(A):
        return spsolve(A=A, b=b)

    elif isinstance(A, np.ndarray) or isinstance(A, np.matrix):
        return np.linalg.solve(a=A, b=b)

    else:
        raise TypeError


def chi2_data(Data, Dobs, CDinv):

    assert sp.issparse(CDinv)
    assert not sp.issparse(Data - Dobs)

    Inan = np.isnan(Data) | np.isnan(Dobs)
    Dobs[Inan] = Data[Inan] = 0.  # so that they do not act on the data misfit

    return np.dot((Data - Dobs), CDinv * (Data - Dobs))  # only because CDinv is known exactly


def sparse_cholesky(A):
    """The input matrix A must be a sparse symmetric positive-definite.
    source https://gist.github.com/omitakahiro/c49e5168d04438c5b20c921b928f1f5d
    """

    n = A.shape[0]
    LU = splu(A, diag_pivot_thresh=0)  # sparse LU decomposition
    plt.figure()
    plt.plot(LU.perm_r)
    plt.figure()
    plt.plot(LU.U.diagonal())
    plt.show()
    if (LU.perm_r == np.arange(n)).all() and (LU.U.diagonal() > 0).all():
        # check the matrix A is positive definite.
        return LU.L.dot(sp.diags(LU.U.diagonal() ** 0.5))
    else:
        raise Exception('The matrix is not positive definite')


def chi2_model(Model, Mprior, CM):

    assert sp.issparse(CM)
    assert not sp.issparse(Model - Mprior)

    cost = np.dot(Model - Mprior, Ainv_dot_b(A=CM, b=Model - Mprior))
    if cost < 0.:
        print('\n\n\n')
        warnings.warn('the model cost was negative (numerical instability?)')
        print('\n\n\n')

        x = Ainv_dot_b(A=CM, b=Model - Mprior)
        plt.figure()
        plt.subplot(211)
        plt.plot(Model - Mprior)
        plt.plot(CM * x)
        plt.subplot(212, sharex=plt.gca())
        plt.plot((Model - Mprior) - CM * x)
        plt.show()


        lu = splu(CM)
        b1 = lu.solve(Model - Mprior)
        b2 = spsolve(CM, Model - Mprior)
        plt.figure()
        plt.plot(b1)
        plt.plot(b2)
        plt.show()


#    L = sparse_cholesky(CM)

    # plt.figure()
    # plt.subplot(211)
    # plt.plot(Model - Mprior)
    # plt.plot(Ainv_dot_b(A=CM, b=Model - Mprior))
    # plt.subplot(212, sharex=plt.gca())
    # plt.plot((Model - Mprior) * Ainv_dot_b(A=CM, b=Model - Mprior))
    # plt.show()
    #

    # use only the diagonal terms of the covariance matrix
    # cost = np.dot(Model - Mprior, Ainv_dot_b(A=sp.diags(CM.diagonal()), b=Model - Mprior))
    return cost


def print_costs(niter, data_cost, model_cost):
    print('iter {:03d}: chi2_data={:<16f} chi2_model={:<16f} chi2={:<16f}'.format(
        niter, data_cost, model_cost, data_cost + model_cost))


def load_last_frechet_derivatives(supertheory, verbose, mapkwargs):
    niter = lastiter()
    while niter >= 0:
        fdfilename = FDFILE.format(niter=niter)
        if os.path.isfile(fdfilename):
            if verbose:
                print('loading {}'.format(fdfilename))
            Gi = sp.load_npz(fdfilename)
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
    #     G0 = sp.load_npz(FDFILE.format(niter=0))
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
        Gi = sp.load_npz(fdfilename)
        return Gi


def rm_nofail(filepath, verbose=True):
    ls = glob.glob(filepath)
    if len(ls):
        for filename in ls:
            if os.path.isfile(filename):
                if verbose:
                    print('removing {}'.format(filename))
                os.remove(filename)
            elif verbose:
                print('could not remove {} (not found)'.format(filename))
    elif verbose:
        print('could not remove {} (no files)'.format(filepath))


def tv23_1(Dobs, Di, CD,
           Mprior, M0, Mi, CM,
           Gi, supertheory,
           CMGiT=None, LUi=None,
           verbose=False):
    """
    Tarantola Valette 1982 eq. 23, modified to control the step size
    :param Dobs:
    :param Di:
    :param CD:
    :param Mprior:
    :param M0:
    :param Mi:
    :param CM:
    :param Gi:
    :param supertheory:
    :return:
    """
    # nan issues
    # nan can occur in the predicted data if a data point is above the cut off period
    # let ignore them
    Inan = np.isnan(Di) | np.isnan(Dobs)
    Dobs[Inan] = Di[Inan] = 0.
    print('computing Dobs - Di + Gi * (Mi - Mprior)...')
    Xi = Dobs - Di + Gi * (Mi - Mprior)
    print('ok')

    if CMGiT is None:
        print('computing CM . Gi.T...')
        CMGiT = CM * Gi.T
        print('ok')

    Ai = None
    if LUi is None:
        print('computing CD + Gi . CM . Gi.T...')
        Ai = CD + Gi * CMGiT
        print('ok')

        print('factorizing CD + Gi . CM . Gi.T...')
        LUi = splu(Ai)
        print('ok')

    print('computing (CD + Gi *. CM . Gi.T)^-1 . (Dobs - Di + Gi . (Mi - Mprior)) ...')
    Aiinv_dot_Xi = LUi.solve(Xi)
    print('ok')

    if Ai is not None:
        error = np.abs(Xi - Ai * Aiinv_dot_Xi).sum()
        print('error on A.x=b : {}'.format(error))

    print('computing CM . Gi.T . (CD + Gi . CM . Gi.T)^-1 . (Dobs - Di + Gi . (Mi - Mprior)) ...')
    KiXi = CMGiT * Aiinv_dot_Xi
    print('ok')

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

    return Dinew, Minew, LUi


def optimize(argv, verbose, mapkwargs):

    if '-h' in argv.keys() or "-help" in argv.keys():
        print(long_help)
        return

    check_keys(argv)

    # === inititate the working directory, clean up previous files, copy node file
    if "-init" in argv.keys():
        nodefilename_in = argv['-init'][0]

        # ========= clear temporary files
        for filename in CLEAR_LIST:
            rm_nofail(filename)

        # ========= get the node file and copy it in the working directory
        if not os.path.isfile(nodefilename_in):
            raise IOError(nodefilename_in)

        nodefile = NodeFile(nodefilename_in)
        nodefile.fill_extraction_files()
        nodefile.write(NODEFILELOCAL)
        assert os.path.isfile(NODEFILELOCAL)

    if "-restart" in argv.keys():
        # ========= clear temporary files except for the initiation files
        m96fileout = M96FILEOUT.format(node="*")
        s96fileout = S96FILEOUT.format(node="*")
        rm_nofail(m96fileout, verbose=verbose)
        rm_nofail(s96fileout, verbose=verbose)

        for niter in range(1, lastiter()+1):
            mfile = MFILE.format(niter=niter)
            dfile = DFILE.format(niter=niter)
            fdfile = FDFILE.format(niter=niter)
            for filename in [mfile, dfile, fdfile]:
                rm_nofail(filename, verbose=verbose)

    # ================================= exit now if possible
    for key in argv.keys():
        if key in ["-prior", "-data", "-fd", "-upd", "-show", "-save"]:
            # continue execution only if one of the options is active
            break
    else:
        import sys
        print('no more options, exiting now')
        sys.exit(0)

    # =================================
    # (re)read the local node file
    nodefile = NodeFileLocal(NODEFILELOCAL)
    nodefile.fill_extraction_files()

    # reload or build the forward problem operator
    n_data_points, n_parameters, \
        superdatacoder, superparameterizer, \
        supertheory = \
            nodefile.build_or_load_forward_problem(
            verbose=verbose, mapkwargs=mapkwargs)

    if "-prior" in argv.keys():
        horizontal_smoothing_distance = argv["-prior"][0]
        vertical_smoothing_distance = argv["-prior"][1]
        trunc_horizontal_smoothing_distance = argv["-prior"][2]
        trunc_vertical_smoothing_distance = argv["-prior"][3]
        lock_half_space = bool(int(argv["-prior"][4]))
        scale_uncertainties = argv["-prior"][5]
        add_uncertainty = argv["-prior"][6]

        # ========= set the prior model and covariance matrix
        if verbose:
            print('get Mprior...')
        Mprior = superparameterizer.get_Mprior()
        if verbose:
            print('done')

            print('get CM_triu ...')
        CM_triu = nodefile.get_CM(
            horizontal_smoothing_distance=horizontal_smoothing_distance,
            vertical_smoothing_distance=vertical_smoothing_distance,
            trunc_horizontal_smoothing_distance=trunc_horizontal_smoothing_distance,
            trunc_vertical_smoothing_distance=trunc_vertical_smoothing_distance,
            lock_half_space=lock_half_space,
            scale_uncertainties=scale_uncertainties,
            add_uncertainty=add_uncertainty,
            norm="L1",
            visual_qc=False,
            mapkwargs=mapkwargs)

        if verbose:
            print('done')

        save_matrix(MPRIORFILE, Mprior, verbose)
        save_matrix(CMFILE, CM_triu, verbose)

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
        if verbose:
            print('loading ...')
        Mprior = np.load(MPRIORFILE)
        M0 = np.load(MFILE.format(niter=0))
        Dobs = np.load(DOBSFILE)

        CM = load_CM(CMFILE)
        CD = sp.load_npz(CDFILE)
        CDinv = sp.load_npz(CDINVFILE)

        LUi = None   # the factorized version of CD + Gi . CM . Gi^T
        if not update_G:
            # find the last version of the frechet derivatives, compute G0 if none is found
            Gi = load_last_frechet_derivatives(supertheory, verbose=verbose, mapkwargs=mapkwargs)
            CMGiT = supertheory.get_CMGiT(CM, Gi, mapkwargs=mapkwargs, verbose=verbose)

        for _ in range(number_of_iterations):
            niter = lastiter()

            if update_G:
                # compute frechet derivatives for this iteration (if not already computed)
                Gi = update_frechet_derivatives(supertheory, verbose=verbose, mapkwargs=mapkwargs)
                CMGiT = supertheory.get_CMGiT(CM, Gi, mapkwargs=mapkwargs, verbose=verbose)
                LUi = None   # force recompute LUi since Gi has been updated

            # ==== save the current model state
            mfilename = MFILE.format(niter=niter)
            dfile = DFILE.format(niter=niter)
            # fdfilename = FDFILE.format(niter=niter)

            Mi = np.load(mfilename)
            Di = np.load(dfile)  # g(Mi)
            # Gi = sp.load_npz(fdfilename)

            # ==== update the current model
            Dinew, Minew, LUi = tv23_1(Dobs, Di, CD,
                       Mprior, M0, Mi, CM,
                       Gi, supertheory,
                       CMGiT=CMGiT,
                       LUi=LUi)

            # ==== save the new model state
            mfilename_new = MFILE.format(niter=niter + 1)
            dfilename_new = DFILE.format(niter=niter + 1)
            # fdfilename_new = FDFILE.format(niter=niter + 1)

            save_matrix(mfilename_new, Minew, verbose)
            save_matrix(dfilename_new, Dinew, verbose)
            # save_matrix(fdfilename_new, Ginew, verbose)

            if verbose:
                data_cost = chi2_data(Data=Dinew, Dobs=Dobs, CDinv=CDinv)
                model_cost = chi2_model(Model=Minew, Mprior=Mprior, CM=CM)
                print_costs(niter+1, data_cost=data_cost, model_cost=model_cost)

    if "-show" in argv.keys():

        modelfiles = np.sort(glob.glob(MFILES))
        datafiles = np.sort(glob.glob(DFILES))

        Dobs = np.load(DOBSFILE)
        Mprior = np.load(MPRIORFILE)
        CM = load_CM(CMFILE)
        CDinv = sp.load_npz(CDINVFILE)
        Mpriorunc = CM.diagonal() ** 0.5

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

            data_cost = chi2_data(Data=Data, Dobs=Dobs, CDinv=CDinv)
            model_cost = chi2_model(Model=Model, Mprior=Mprior, CM=CM)

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

