import numpy as np
import matplotlib.pyplot as plt


"""
goal : find the model whose power of 3 fits some data
using linearized approach with a non diagonal covariance matrix to smooth the model
"""


def stairs(x, y):
    dx = np.median(x[1:] - x[:-1])
    x_ = np.hstack((x[0] - dx / 2,
                    np.repeat(0.5 * (x[1:] + x[:-1]), 2),
                    x[-1] + dx / 2.))
    y_ = np.repeat(y, 2)
    return x_, y_

def lininvsmooth(xi, di, si, xmodel, modelprior, sigmamodel,
              correlation_length=0.1):

    n = len(di)
    assert len(xi) == len(si) == n

    m = len(modelprior)
    assert len(xmodel) == m

    Cd = np.diag(si ** 2.0)
    if correlation_length == 0.:
        # no smoothing
        Cm = np.diag(sigmamodel ** 2.0)
    elif False:
        # use a exponential smoothing
        Cm = sigmamodel[:, np.newaxis] * sigmamodel * np.exp(
            -np.abs((xmodel[:, np.newaxis] - xmodel) / correlation_length))
    else:
        # use a gaussian smoothing
        Cm = sigmamodel[:, np.newaxis] * sigmamodel * \
             np.exp(-0.5 * ((xmodel[:, np.newaxis] - xmodel) / correlation_length) ** 2.)

    # xmodel_filt = 0.5 * (xmodel[1:] + xmodel[:-1])
    xmodel_filt = xmodel

    def filt(model):

        # return ((model[1:] - model[:-1]) / (xmodel[1:] - xmodel[:-1]))
        return model ** 3.0

    def g(model):
        # compute the derivative of model
        model_filt = filt(model)

        # find the values at xi
        dcalc = np.interp(xi, xp=xmodel_filt, fp=model_filt)

        return dcalc

    def FD(model):
        fd = np.zeros((n, m), float)
        gm = g(model)
        for j in range(m):
            model_j = model.copy()
            model_j[j] += 0.01

            fd[:, j] = (g(model_j) - gm) / 0.01

        return fd

    modeli = modelprior
    models = [modeli.copy()]

    Cdinv = np.linalg.inv(Cd)  # np.diag(si ** -2.0)
    Cminv = np.linalg.inv(Cm)
    for _ in range(100):

        Gi = FD(modeli)
        if n <= m:
            # over determined problem
            CmGiT = np.dot(Cm, Gi.T)
            Siinv = np.linalg.inv(Cd + np.dot(Gi, CmGiT))
            # Siinv = np.diag(np.diag(Siinv))
            Hi = np.dot(CmGiT, Siinv)

        else:
            # under determined problem
            GiTCdinv = np.dot(Gi.T, Cdinv)
            Siinv = np.linalg.inv(np.dot(GiTCdinv, Gi) + Cminv)
            Hi = np.dot(Siinv, GiTCdinv)

        Xi = di - g(modeli) + np.dot(Gi, (modeli - models[0]))
        modeli = models[0] + np.dot(Hi, Xi)

        models.append(modeli.copy())
        if len(models) > 1:
            if ((models[-1] - models[-2]) ** 2.0 / float(m)).sum() ** 0.5 < 0.01:
                # convergence reached
                break

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
        x__, y__ = stairs(xmodel_filt, filt(models[i]))
        # plt.plot(xmodel_filt, filt(models[i]), 'bo', alpha=0.2)
        plt.plot(x__, y__, 'b-', alpha=alpha, linewidth=linewidth, label=label_fit)


if __name__ == '__main__':
    for test in [1, 2]:
        if test == 1:
            # over determined problem
            xi = np.array([0., 1., 2., 3.])
            di = np.array([-2., 1.5, 2.5, -3.5])
            si = np.array([0.1, 0.1, 0.1, 0.1])
            m = 100
        elif test == 2:
            # under determined problem
            xi = np.sort(np.random.rand(100)) * 3.
            di = np.cos(5 * xi) + 4 * xi - 2.5 + 0.5 * np.random.randn(len(xi))
            si = 0.5 * np.ones_like(xi)
            m = 10
        xmodel = np.linspace(xi[0]-0.05, xi[-1]+0.05, m)
        modelprior = np.linspace(-0.5, 0.5, m) # + 0.1 * np.random.randn(m)
        sigmamodel = 0.5 * np.ones(m)
        correlation_length = 0.3

        plt.figure()
        lininvsmooth(xi, di, si, xmodel, modelprior, sigmamodel, correlation_length)
        plt.gca().set_ylim((min(di)-1, max(di)+1))
        plt.legend()
    plt.ion()
    plt.show()
    input('pause')
