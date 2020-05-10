from __future__ import print_function
from builtins import input
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import issparse, diags, csr_matrix, csc_matrix
from scipy.sparse.linalg import spsolve
# from scipy.sparse.linalg import inv as sparse_inv
# from numpy.linalg import inv


def tv23(mprior, mn, CM,
        dobs, dn, CD,
        Gn, mun=1.0):

    n_param = len(mprior)
    n_data = len(dobs)

    # seems stable for both under and overdetermined problems

    assert len(mn) == n_param
    assert CM.shape == (n_param, n_param)
    assert len(dn) == len(dobs)
    assert CD.shape == (n_data, n_data)
    assert Gn.shape == (n_data, n_param)

    # === operators that must be a np.matrix or a scipy.spmatrix
    if not issparse(Gn):
        Gn = np.asmatrix(Gn)
    if not issparse(CM):
        CM = np.asmatrix(CM)
    if not issparse(CD):
        CD = np.asmatrix(CD)

    # ===
    for mat in [mprior, mn, dobs, dn]:
        assert not issparse(mat)
        assert isinstance(mat, np.ndarray)  # or isinstance(mat, np.matrix)
        assert mat.ndim == 1
    mprior = np.asmatrix(mprior[:, np.newaxis])
    mn = np.asmatrix(mn[:, np.newaxis])
    dobs = np.asmatrix(dobs[:, np.newaxis])
    dn = np.asmatrix(dn[:, np.newaxis])

    # === Kn
    CMGT = CM * Gn.T
    Sn = CD + Gn * CM * Gn.T

    if 0:
        # use only the diag terms of Sn
        assert isinstance(Sn, np.matrix)
        Sn = diags(np.diag(np.asarray(Sn)), shape=(n_data, n_data))
        # WRONG : DO NOT WORK

    # if issparse(Sn):
    #     Sninv = sparse_inv(Sn)
    # elif isinstance(Sn, np.matrix):
    #     Sninv = inv(Sn)
    # else:
    #     raise TypeError
    #
    # if True:
    #     import warnings
    #     warnings.warn('test')
    #     # use only the diag terms of Sn
    #     assert isinstance(Sninv, np.matrix)
    #     Sninv = diags(np.diag(np.asarray(Sninv)), shape=(n_data, n_data))
    #     # WEIRD works only for n_data < n_param

    assert issparse(Sn) or isinstance(Sn, np.matrix)
    #Kn = CMGT * Sninv
    #dm = Kn * (dobs - dn + Gn * (mn - mprior))
    #dm = CMGT * Sninv * (dobs - dn + Gn * (mn - mprior))

    dm = CMGT * Ainv_dot_b(A=Sn, b=(dobs - dn + Gn * (mn - mprior)))[:, np.newaxis]

    assert issparse(dm) or isinstance(dm, np.matrix)

    mnpp = mn + mun * (mprior + dm - mn)

    return np.asarray(mnpp).flat[:]


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
        Ainvb = spsolve(A=A, b=b)

    elif isinstance(A, np.ndarray) or isinstance(A, np.matrix):
        Ainvb = np.linalg.solve(a=A, b=b)

    else:
        raise TypeError

    return np.asarray(Ainvb).flat[:]


if __name__ == '__main__':
    if 1:
        # over determined problem
        # plus je met de la covariance plus l'egalite A20 de TV est fausse
        xobs = np.array([0., 1., 2., 3.])
        dobs = np.array([-2., 1.5, 2.5, -3.5])
        dunc = np.array([0.1, 0.1, 0.1, 0.1])
        m = 100
    else:
        # under determined problem
        xobs = np.sort(np.random.rand(100)) * 3.
        dobs = np.cos(5 * xobs) + 4 * xobs - 2.5 + 0.5 * np.random.randn(len(xobs))
        dunc = 0.5 * np.ones_like(xobs)
        m = 10

    n = len(dobs)
    CD = diags(dunc ** 2.0, format="csc", shape=(n, n))
    CDinv = diags(dunc ** -2.0, format="csc", shape=(n, n))

    xmodel = np.linspace(xobs[0] - 0.05, xobs[-1] + 0.05, m)
    modelprior = np.linspace(-0.5, 0.5, m)  # + 0.1 * np.random.randn(m)
    sigmamodel = 0.5 * np.ones(m)
    correlation_length = 0.3

    # ================ model covariance
    if correlation_length == 0.:
        # no smoothing
        CM = np.diag(sigmamodel ** 2.0)
    elif True:
        # use a exponential smoothing
        CM = sigmamodel[:, np.newaxis] * sigmamodel * np.exp(
            -np.abs((xmodel[:, np.newaxis] - xmodel) / correlation_length))
    else:
        # use a gaussian smoothing
        CM = sigmamodel[:, np.newaxis] * sigmamodel * \
             np.exp(-0.5 * ((xmodel[:, np.newaxis] - xmodel) / correlation_length) ** 2.)

    # CMinv = np.linalg.inv(CM)

    # ================
    def filt(model):

        # return ((model[1:] - model[:-1]) / (xmodel[1:] - xmodel[:-1]))
        return model ** 3.0

    def g(model):
        # compute the derivative of model
        model_filt = filt(model)

        # find the values at xi
        dcalc = np.interp(xobs, xp=xmodel, fp=model_filt)

        return dcalc

    def FD(model):
        fd = np.zeros((n, m), float)
        gm = g(model)
        for j in range(m):
            model_j = model.copy()
            model_j[j] += 0.01

            fd[:, j] = (g(model_j) - gm) / 0.01

        return fd

    # ================
    mn = modelprior
    dn = g(mn)
    Gn = FD(mn)
    # raise ValueError(mn.shape, dn.shape, Gn.shape)

    models = [mn]
    for _ in range(100):
        mn = tv23(mprior=modelprior,
             mn=mn,
             CM=CM,
             # CMinv=CMinv,
             dobs=dobs,
             dn=dn,
             CD=CD,
             # CDinv=CDinv,
             Gn=Gn,
             mun=0.6)
        dn = g(mn)
        Gn = FD(mn)
        models.append(mn)

    def stairs(x, y):
        dx = np.median(x[1:] - x[:-1])
        x_ = np.hstack((x[0] - dx / 2,
                        np.repeat(0.5 * (x[1:] + x[:-1]), 2),
                        x[-1] + dx / 2.))
        y_ = np.repeat(y, 2)
        return x_, y_

    xi = xobs
    di = dobs
    si = dunc
    plt.plot(xi, di, 'ko', label="data")
    for i in range(n):
        plt.plot([xi[i],xi[i]], [di[i]-si[i], di[i]+si[i]], "k_-")

    for i in range(len(models)):
        linewidth = 1
        alpha = 0.2
        label_model = None
        label_fit = None
        if i == len(models) - 1:
            linewidth = 3
            alpha = 1
            label_model = "solution"
            label_fit = "data fit"
        elif i == 0:
            linewidth = 1
            alpha = 1
            label_model = "prior"

        x_, y_ = stairs(xmodel, models[i])
        # plt.plot(xmodel, models[i], 'ro', alpha=0.2)
        plt.plot(x_, y_, 'r-', alpha=alpha, linewidth=linewidth, label=label_model)
        x__, y__ = stairs(xmodel, filt(models[i]))
        # plt.plot(xmodel_filt, filt(models[i]), 'bo', alpha=0.2)
        plt.plot(x__, y__, 'b-', alpha=alpha, linewidth=linewidth, label=label_fit)

    plt.legend()
    plt.ion()
    plt.show()
    input('pause')