import warnings
import numpy as np
import matplotlib.pyplot as plt
from srfpython.standalone.stdout import waitbarpipe
from srfpython.standalone.multipro8 import Job, MapAsync
from srfpython.HerrMet.theory import Theory
from srfpython.depthdisp.dispcurves import surf96reader_from_arrays
from srfpython.HerrMet.paramfile import load_paramfile
from srfpython.HerrMet.datacoders import Datacoder_log, makedatacoder
from srfpython.standalone.multipro8 import Job, MapSync, MapAsync
from srfpython.standalone.stdout import waitbarpipe
from srfpython.utils import Timer
import scipy.sparse as sp
from scipy.sparse import linalg as splinalg


# ========== g
class ForwardOperator(object):

    def __init__(self,
                 parameterizer_strings, datacoder_strings,
                 nx, ny, nz,
                 verbose=True, **mapkwargs):
        self.mapkwargs = mapkwargs
        self.nx, self.ny, self.nz = nx, ny, nz

        def job_generator():
            ls = zip(parameterizer_strings, datacoder_strings)
            for nnode, (ps, ds) in enumerate(ls):
                yield Job(nnode, parameterizer_string=ps, datacoder_string=ds)

        def job_handler(nnode, parameterizer_string, datacoder_string):
            parameterizer = load_paramfile(parameterizer_string, verbose=False)[0]
            datacoder = makedatacoder(datacoder_string, which=Datacoder_log)

            theory = Theory(parameterizer=parameterizer, datacoder=datacoder)

            return nnode, theory

        wb = None
        if verbose:
            wb = waitbarpipe('build g')

        theorys = []
        with MapSync(job_handler, job_generator(), **self.mapkwargs) as ma:
            for jobid, (nnode, theory), _, _ in ma:
                theorys.append(theory)

                if verbose:
                    wb.refresh(jobid / float(len(parameterizer_strings)))

        if verbose:
            wb.close()

        self.theorys = np.array(theorys, dtype=object).reshape((ny, nx))

    def __call__(self, M, verbose=True):
        M = M.reshape((nz, ny, nx))
        theorys = self.theorys

        def job_generator():
            for iy in range(ny):
                for jx in range(nx):  # order matters!!!!
                    nnode = iy * nx + jx
                    yield Job(nnode, theorys[iy, jx], M[:, iy, jx])

        def job_handler(nnode, theory, m):
            data = theory(m=m)
            return nnode, data

        wb = None
        if verbose:
            wb = waitbarpipe('g(m)')

        Data = np.zeros_like(Dobs)
        with MapAsync(job_handler, job_generator(), **self.mapkwargs) as ma:

            for nnode, (nnode, data), _, _ in ma:
                ib = Nobs[:nnode].sum()
                ie = ib + Nobs[nnode]

                Data[ib: ie] = data

                if verbose:
                    wb.refresh(nnode / float(nx * ny))

        if verbose:
            wb.close()

        return Data  # warning : Data means encoded data

    def frechet_derivatives(self, M, verbose=True):

        M = M.reshape((nz, ny, nx))
        theorys = self.theorys

        def job_generator():
            for iy in range(ny):
                for jx in range(nx): # order matters!!!!
                    nnode = iy * nx + jx
                    yield Job(nnode, theorys[iy, jx], M[:, iy, jx])

        def job_handler(nnode, theory, m):
            fd = theory.frechet_derivatives(m=m)
            return nnode, fd

        wb = None
        if verbose:
            wb = waitbarpipe('dg/dm(m)')

        # FDs = []
        rows = []
        cols = []
        dats = []
        with MapSync(job_handler, job_generator(), **self.mapkwargs) as ma:  # order matters
            for jobid, (nnode, fd), _, _ in ma:
                # FDs.append(fd)
                iy = nnode // nx  # line number of the node in the surface plane = "latitude number"
                ix = nnode % nx  # column number of the node in the surface plane = "longitude number"

                # _cols, _rows = np.meshgrid(np.arange(nz), np.arange(nper))
                # cols += list(nnode * nz + _cols.flat[:])
                # rows += list(nnode * nper + _rows.flat[:])

                # iz = depth number
                iz, ida = np.meshgrid(np.arange(nz), np.arange(Nobs[nnode]))
                # convert indexs into positions in the model (_cols) ad data (_rows) spaces
                _cols = iz.flat[:] * (nx * ny) + iy * nx + ix
                _rows = list(Nobs[:nnode].sum() + ida.flat[:])
                rows += list(_rows)
                cols += list(_cols)
                dats += list(fd.flat[:])

                if verbose:
                    wb.refresh(nnode / float(nx * ny))

        if verbose:
            wb.close()

        G = sp.csc_matrix((dats, (rows, cols)), shape=(len(Dobs), len(Mprior)), dtype=float)
        # plt.figure()
        # plt.imshow(G.A)
        # plt.show()
        return G


# ========== CM
class ModelSmoother(object):
    def __init__(self,
                 x, y, zmid,
                 xflat, yflat, zmidflat,
                 Lh, Lv):
        self.Lh = Lh
        self.Lv = Lv
        self.x, self.y, self.zmid = x, y, zmid
        self.xflat = xflat
        self.yflat = yflat
        self.zmidflat = zmidflat

        self.nx, self.ny, self.nz = len(x), len(y), len(zmid)
        self.shape = (self.nx * self.ny * self.nz, self.nx * self.ny * self.nz)

    def row(self, i, columns=slice(None, None, None)):
        Lh, Lv = self.Lh, self.Lv
        xflat, yflat, zflat = self.xflat, self.yflat, self.zmidflat

        d2norm = ((xflat[columns] - xflat[i]) / Lh) ** 2.0 \
               + ((yflat[columns] - yflat[i]) / Lh) ** 2.0 \
               + ((zflat[columns] - zflat[i]) / Lv) ** 2.0
        s_line = np.exp(-0.5 * d2norm)
        return s_line / s_line.sum()

    def col(self, j, rows=slice(None, None, None)):
        # CM.T = CM
        return self.row(i=j, columns=rows)

    def sparse_row(self, i, trunc=2.):
        x, y, z = self.x, self.y, self.zmid
        nx, ny, nz = self.nx, self.ny, self.nz
        Lh, Lv = self.Lh, self.Lv
        xflat, yflat, zflat = self.xflat, self.yflat, self.zmidflat

        ix = np.abs(x - xflat[i]) < trunc * Lh
        iy = np.abs(y - yflat[i]) < trunc * Lh
        iz = np.abs(zmid - zflat[i]) < trunc * Lv

        I = np.arange(nz * ny * nx).reshape((nz, ny, nx))[iz, ...][:, iy, :][..., ix].flat[:]

        cols = np.arange(self.shape[1])[I]
        rows = i * np.ones_like(cols)

        d2norm = (np.abs(xflat[I] - xflat[i]) / Lh) ** 2.0 \
                +(np.abs(yflat[I] - yflat[i]) / Lh) ** 2.0 \
                +(np.abs(zflat[I] - zflat[i]) / Lv) ** 2.0
        vals = np.exp(-0.5 * d2norm)
        vals /= vals.sum()

        return rows, cols, vals

    def dot(self, b, trunc=0.):

        assert b.ndim == 1
        assert len(b) == self.shape[1]
        assert isinstance(b, np.ndarray)
        self_dot_b = np.zeros_like(b)

        def job_generator():
            for i in range(len(b)):
                yield Job(i)

        def job_handler(i):
            if trunc <= 0.:
                return i, (self.row(i) * b).sum()
            else:
                rows, cols, vals = self.sparse_row(i, trunc=trunc)
                return i, (vals * b[cols]).sum()

        wb = waitbarpipe('dot product')
        with MapAsync(job_handler, job_generator()) as ma:
            for _, (i, v), _, _ in ma:
                self_dot_b[i] = v
                wb.refresh(i / float(self.shape[0]))
            wb.close()

        return self_dot_b

    @property
    def A(self):
        self_dense = np.zeros(self.shape, float)
        for i in range(self.shape[0]):
            self_dense[i, :] = self.row(i)
        return self_dense


class ModelCovarianceMatrix(ModelSmoother):
    # shape = (len(Mprior), len(Mprior))

    def diagonal(self, offset=0):
        if offset == 0:
            return Munc ** 2.0
        elif 0 < offset < self.shape[1]:
            raise NotImplementedError
        elif offset < 0:
            return self.diagonal(offset=-offset)
        else:
            raise ValueError(offset)

    def row(self, i, columns=slice(None, None, None)):
        Lh, Lv = self.Lh, self.Lv
        xflat, yflat, zflat = self.xflat, self.yflat, self.zmidflat

        d2norm = ((xflat[columns] - xflat[i]) / Lh) ** 2.0 \
                +((yflat[columns] - yflat[i]) / Lh) ** 2.0 \
                +((zflat[columns] - zflat[i]) / Lv) ** 2.0

        return Munc[columns] * Munc[i] * np.exp(-0.5 * d2norm)

    def sparse_row(self, i, trunc=2.):
        x, y, z = self.x, self.y, self.zmid
        nx, ny, nz = self.nx, self.ny, self.nz
        Lh, Lv = self.Lh, self.Lv
        xflat, yflat, zflat = self.xflat, self.yflat, self.zmidflat

        ix = np.abs(x - xflat[i]) < trunc * Lh
        iy = np.abs(y - yflat[i]) < trunc * Lh
        iz = np.abs(zmid - zflat[i]) < trunc * Lv

        I = np.arange(nz * ny * nx).reshape((nz, ny, nx))[iz, ...][:, iy, :][..., ix].flat[:]

        cols = np.arange(self.shape[1])[I]
        rows = i * np.ones_like(cols)

        d2norm = (np.abs(xflat[I] - xflat[i]) / Lh) ** 2.0 \
                +(np.abs(yflat[I] - yflat[i]) / Lh) ** 2.0 \
                +(np.abs(zflat[I] - zflat[i]) / Lv) ** 2.0
        vals = Munc[I] * Munc[i] * np.exp(-0.5 * d2norm)

        return rows, cols, vals

    def sparse(self, trunc=2.0):

        def job_generator():
            for i in range(self.shape[0]):
                yield Job(i)

        def job_handler(i):
            _rows, _cols, _vals = self.sparse_row(i, trunc=trunc)
            return _rows, _cols, _vals

        rows = []
        cols = []
        vals = []
        with MapAsync(job_handler, job_generator(), Nworkers=40) as ma:
            for _, (_rows, _cols, _vals), _, _ in ma:
                rows.append(_rows)
                cols.append(_cols)
                vals.append(_vals)

        rows = np.hstack(rows)
        cols = np.hstack(cols)
        vals = np.hstack(vals)

        return sp.csc_matrix((vals, (rows, cols)))


if __name__ == '__main__':

    # ========== Load inputs
    with np.load('mprior.npz') as loader:
        x = loader['x']
        y = loader['y']
        zmid = loader['zmid']
        xflat = loader['xflat']
        yflat = loader['yflat']
        zflat = zmidflat = loader['zmidflat']
        Lh = loader['Lh']
        Lv = loader['Lv']
        Mprior = loader['Mprior']
        Munc = loader['Munc']
        damping = loader['damping']
        parameterizer_strings = loader['parameterizer_strings']

    with np.load('dobs.npz') as loader:
        Nobs = loader['Nobs']
        Waves = loader['Waves']
        Types = loader['Types']
        Modes = loader['Modes']
        Freqs = loader['Freqs']
        Dobs = loader['Dobs']
        Dunc = loader['Dunc']
        datacoder_strings = loader['datacoder_strings']

    nx, ny, nz = len(x), len(y), len(zmid)
    nobs = len(Dobs)

    Munc *= np.sqrt(float(len(Mprior)) / float(len(Dobs))) / damping  # equilibrate the costs, apply damping here

    CD = sp.diags(Dunc ** 2., format="csc", shape=(len(Dobs), len(Dobs)))
    CDinvdiag = Dunc ** -2.

    smoother = ModelSmoother(
        x=x, y=y, zmid=zmid,
        xflat=xflat, yflat=yflat, zmidflat=zflat,
        Lh=Lh, Lv=Lv)
    CM = ModelCovarianceMatrix(
        x=x, y=y, zmid=zmid,
        xflat=xflat, yflat=yflat, zmidflat=zflat,
        Lh=Lh, Lv=Lv)
    g = ForwardOperator(
        parameterizer_strings=parameterizer_strings,
        datacoder_strings=datacoder_strings,
        nx=nx, ny=ny, nz=nz,
        Nworkers=8)

    Dprior = g(Mprior)
    np.save('Dprior.npy', Dprior)  # before smoothing

    plt.figure()
    plt.plot(Dobs)
    plt.plot(Dprior)
    plt.show()

    M0 = smoother.dot(Mprior, trunc=4.)

    D0 = g(M0)
    np.save('M0.npy', M0)
    np.save('D0.npy', D0)

    MI = M0.copy()
    DI = D0.copy()

    with Timer('CM_sparse'):
        CM_sparse = CM.sparse(trunc=4.)
    with Timer('CM_lu'):
        CM_sparse_lu = splinalg.splu(CM_sparse)

    # plt.figure()
    # plt.imshow(CM_sparse.A)
    # plt.show()
    # exit()

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

        np.save('data_costs.npy', np.asarray(data_costs))
        np.save('model_costs.npy', np.asarray(model_costs))
        np.save('total_costs.npy', np.asarray(model_costs) + np.asarray(data_costs))
        np.save('Msol.npy', MI)
        np.save('Dsol.npy', DI)

        if np.abs(DM).max() <= 0.01:
            print('convergence achieved')
            break