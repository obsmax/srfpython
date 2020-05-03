import matplotlib.pyplot as plt
import numpy as np

"""
tester for linear inversion with smoothing
"""


def smoothinv(xi, di, si,
              xmodel, modelprior, sigmamodel,
              correlation_length=0.1):
    m = len(xmodel)
    n = len(xi)
    assert len(modelprior) == len(sigmamodel) == m

    G = np.zeros((n, m), float)
    for i in range(n):
        j = np.argmin(abs(xmodel - xi[i]))
        G[i, j] = 1.0

    Cd = np.diag(si ** 2.0)
    if False:
        Cm = np.diag(sigmamodel ** 2.0)
    elif False:
        Cm = sigmamodel[:, np.newaxis] * sigmamodel * np.exp(
            -np.abs((xmodel[:, np.newaxis] - xmodel) / correlation_length))
    else:
        Cm = sigmamodel[:, np.newaxis] * sigmamodel * np.exp(-0.5 * ((xmodel[:, np.newaxis] - xmodel) / correlation_length) ** 2.)

    CmGT = np.dot(Cm, G.T)
    Sinv = np.linalg.inv(Cd + np.dot(G, CmGT))

    model = modelprior + np.dot(CmGT, np.dot(Sinv, (di - np.dot(G, modelprior))))

    plt.plot(xi, di, 'ko', label="data")
    for i in range(n):
        plt.plot([xi[i],xi[i]], [di[i]-si[i], di[i]+si[i]], "k_-")

    xmodel_ = np.hstack((xmodel[0], np.repeat(xmodel[1:], 2), xmodel[-1] + xmodel[1] - xmodel[0]))

    plt.plot(xmodel_, np.repeat(modelprior, 2), 'r-', label="prior")
    plt.fill_between(xmodel_,
                     np.repeat(modelprior - sigmamodel, 2),
                     np.repeat(modelprior + sigmamodel, 2),
                     color='r', alpha=0.1)
    plt.plot(xmodel_, np.repeat(model, 2), 'b-', label="solution")
    return model


if __name__ == '__main__':
    xi = np.array([0., 1., 2., 3.])
    di = np.array([-1., 1., 1., -1.])
    si = np.array([0.1, 0.1, 0.1, 0.1])
    m = 100
    xmodel = np.linspace(xi[0], xi[-1], m)
    modelprior = np.linspace(-0.5, 0.5, m)
    sigmamodel = 0.25 * np.ones(m)

    correlation_length = 0.1
    smoothinv(xi, di, si, xmodel, modelprior, sigmamodel, correlation_length)

    correlation_length = 2.0
    smoothinv(xi, di, si, xmodel, modelprior, sigmamodel, correlation_length)
    # modelprior = smoothinv(xi, di, si, xmodel, modelprior, sigmamodel, correlation_length)
    # modelprior = smoothinv(xi, di, si, xmodel, modelprior, sigmamodel, correlation_length)
    # modelprior = smoothinv(xi, di, si, xmodel, modelprior, sigmamodel, correlation_length)
    # modelprior = smoothinv(xi, di, si, xmodel, modelprior, sigmamodel, correlation_length)
    # modelprior = smoothinv(xi, di, si, xmodel, modelprior, sigmamodel, correlation_length)


    # sigmamodel /= 2.
    # smoothinv(xi, di, si, xmodel, modelprior, sigmamodel, correlation_length)

    plt.legend()
    plt.show()