from tetedenoeud.inversion.tarantola import inverse
import numpy as np

"""
inversion linearisee iterative
cf cahier 14, 25/02/2017

method1
misfit   = (di - diobs)
model    = (mj - mjapr)
gradient = finite difference with step eps 

method2 :
not implemented

method3
misfit = (di - diobs) / (abs(di) + eps)
model  = (mj - mjapr) / (abs(mj) + eps)
gradient = finite difference with step pertu * abs(mj[j]) + eps

"""

def minmax(X):
    if hasattr(X, "min"):
        return X.min(), X.max()
    return min(X), max(X)


def l2a(*args):
    return [np.asarray(_) for _ in args]

# -----------------------------------
def K(g, m, gm=None, locked=None, pertu=0.05, eps=1.e-6, method=1):
    # numerical gradient / sensitivity kernel
    # assert np.all(m != 0.)

    if gm is None: gm = g(m)
    if locked is None:
        locked = np.zeros(len(m), bool)
    else:
        assert len(locked) == len(m)
        assert not np.all(locked)
    # assert np.all(gm != 0.)

    Kout = np.zeros((len(gm), len(m)), float)
    for j in xrange(len(m)):
        if locked[j]:
            # jth parameter of m is locked, keep derivative at 0
            continue
        if method == 1:
            dmj = np.zeros_like(m)
            dmj[j] = eps
        elif method == 3:
            dmj = np.zeros_like(m)
            dmj[j] = abs(m[j]) * pertu + eps
        else:
            raise Exception('')

        gmj = g(m + dmj)

        #        if method == 1:
        #            #perturbation has no effect, cannot compute gradient
        #            while ((gmj - gm) == 0.).all():
        #                dmj[j] *= 1.01
        #                gmj = g(m + dmj)

        if method == 1:
            Kout[:, j] = (gmj - gm) / dmj[j]
        elif method == 3:
            Kout[:, j] = (gmj - gm) / (abs(gm) + eps)
            Kout[:, j] *= (abs(m[j]) + eps) / dmj[j]
        else:
            raise Exception('')
    return Kout


# -----------------------------------
class LinInv(object):
    def __init__(self, g, d, CDinv, mapr, CMinv, damping=1.0, locked=None, pertu=0.05, eps=1.e-6, method=1):
        """
        :param g: a callable that gets a data array (d) and returns a model array m
        :param d: the target data array
        :param CDinv: inverted covariance matrix(2D) on the data,
            or inverse of the diagonal terms of the cov matrix (1D),
            or inverse the covariance (scalar)
        :param mapr: apriori model array
        :param CMinv: inverted covariance matrix(2D) on the model,
            or inverse of the diagonal terms of the cov matrix (1D),
            or inverse the covariance (scalar)
        :param damping: damping coefficient
        :param locked: boolean array, use it to prevent prameters from beiing inverted
        :param pertu: relative amplitude of the perturbation
        :param eps: small number to avoid division by 0
        :param method: 1 = we invert the difference between d and g(m)
                       3 = we invert the relative difference (d - g(m)) / (g(m) + eps)
        """
        self.g = g
        self.d = d
        self.mapr = mapr
        self.CDinv = CDinv  # diagonal terms (1D array)
        self.CMinv = CMinv  # diagonal terms (1D array) or full covariance matrix (2D array)
        self.pertu = pertu
        self.eps = eps
        self.method = method
        self.damping = damping
        if locked is None:
            self.locked = np.zeros(len(mapr), bool)
        else:
            assert len(locked) == len(mapr)
            self.locked = locked
        # -------------------------------
        self.mt = mapr
        self.dt = None
        self.G = None
        self.mb = self.mt
        self.db = None
        self.state = "init"

    # ----------------------------------
    def lock(self, j):  # lock jth parameter (integer or array or slice)
        self.locked[j] = True

    def unlock(self, j):  # unlock jth parameter (integer or array or slice)
        self.locked[j] = False

    # -----------------------------------
    def grad(self):
        if self.state != "grad":
            assert self.state == "upd"
            # assert np.all(self.mt != 0.)
            # assert np.all(self.dt != 0.)
            self.G = K(self.g, m=self.mt, gm=self.dt, locked=self.locked, pertu=self.pertu, eps=self.eps,
                       method=self.method)
            # print "G", self.G
            self.state = "grad"
        return self.mt, self.dt, self.G

    # -----------------------------------
    def inv(self):
        if self.state != "inv":
            assert self.state == "grad"

            if self.method == 1:
                D = (self.d - self.dt)
                Mprior = np.zeros_like(self.mt)
                Igood = (~np.isnan(D)) & (~np.isinf(D))
                assert Igood.any()
                # ------
                D = D[Igood]
                CDinv = self.CDinv[Igood]
                CMpriorinv = self.CMinv
                G = self.G[Igood, :]
            elif self.method == 3:
                D = (self.d - self.dt) / (abs(self.dt) + self.eps)
                Igood = (~np.isnan(D)) & (~np.isinf(D))
                if (~Igood).any(): raise NotImplementedError('nan data not implemeted for method 3')
                CDinv = self.CDinv[Igood] * (
                        abs(self.dt[Igood]) + self.eps) ** +2.  # self.CD = 1D array for diagonal terms
                Mprior = np.zeros_like(self.mt)
                CMpriorinv = self.CMinv * (
                        abs(self.mt) + self.eps) ** +2.  # diagonal of the covariance array (std ** 2.)
                G = self.G[Igood, :]
            else:
                raise Exception('')
            # print "CMprior", CMprior
            Mpost, CMpost, R, Dpost, CDpost, Dprior, CDprior = inverse(D, CDinv, G, Mprior,
                                                                       CMpriorinv * self.damping ** +2.)

            if self.method == 1:
                self.mb = self.mt + Mpost
                self.db = np.zeros_like(self.dt) * np.nan
                self.db[Igood] = self.dt[Igood] + Dpost

            elif self.method == 3:
                self.mb = self.mt + (abs(self.mt) + self.eps) * Mpost  # self.mt * (Mpost + 1.)
                self.db = self.dt + (abs(self.dt) + self.eps) * Dpost  # self.dt * (Dpost + 1.)

            else:
                raise Exception('')
            self.state = "inv"

        return self.mb, self.db  # inverted model and data, reduced data misfit on next iteration

    # -----------------------------------
    def upd(self):
        "move current model to the inverted one"
        if self.state != "upd":
            assert self.state in ["init", "inv"]
            self.mt = self.mb
            self.dt = self.g(self.mt)  # predicited data

            # ---------------------
            if self.method == 1:
                D = (self.d - self.dt)
                Igood = (~np.isnan(D)) & (~np.isinf(D))
                self.chi2red = np.sum(D[Igood] ** 2. * self.CDinv[Igood]) / float(Igood.sum())
            elif self.method == 3:
                D = (self.d - self.dt) / (abs(self.dt) + self.eps)
                Igood = (~np.isnan(D)) & (~np.isinf(D))
                if (~Igood).any(): raise NotImplementedError('method 3 not implemented for nan data')
                CDinv = self.CDinv * (abs(self.dt) + self.eps) ** +2.  # self.CD = 1D array for diagonal terms
                D = (self.d - self.dt) / (abs(self.dt) + self.eps)
                self.chi2red = np.sum((D) ** 2. * CDinv) / float(len(D))
            # ---------------------

            self.state = "upd"
        return self.mt, self.dt, self.chi2red


# -----------------------------------
def demo2D():
    # theory
    def g(m):
        """user defined"""
        x, y = m
        x0, y0 = 2., 2.
        if x == x0 and y == y0:
            d0 = 1.0
        else:
            DDD = (((x - x0) / 1.0) ** 2. + ((y - y0) / 1.0) ** 2.) ** 0.5
            d0 = DDD  # np.sin(DDD) / DDD
        d = np.asarray([d0])
        return d

    # -----------------------------------
    plt.figure(figsize=(12, 6))
    ax0 = gcf().add_subplot(121)

    x = np.linspace(-6., 6., 100)
    y = np.linspace(-6., 6., 110)
    X, Y = np.meshgrid(x, y)
    Z = np.asarray([g(np.asarray([xx, yy]))[0] for xx, yy in zip(X.flat[:], Y.flat[:])]).reshape(X.shape)
    gca().pcolormesh(x, y, Z)
    gca().set_aspect(1.)

    # -----------------------------------
    # data
    d = np.array([7.5])
    CDinv = np.array([1.0]) ** -2.
    mapr = np.array([1., 0.0])
    CMinv = np.array([0.5, 0.5]) ** -2.

    # -----------------------------------
    ax1 = gcf().add_subplot(122)
    ax0.contour(X, Y, Z, levels=[d[0]], colors="k")
    ax0.contour(X, Y, Z, levels=[d[0] - CDinv[0] ** -.5, d[0] + CDinv[0] ** -.5], colors="k", linestyles="--")
    ax0.plot(mapr[0], mapr[1], 'ro')

    L = LinInv(g=g, d=d, CDinv=CDinv, mapr=mapr, CMinv=CMinv)
    mt, dt, chi2red = L.upd()
    for n in xrange(10):
        mt, dt, G = L.grad()
        mb, db = L.inv()
        mt, dt, chi2red = L.upd()
        ax0.plot(L.mt[0], L.mt[1], 'w*')
        ax1.plot(n, L.chi2red, 'ko')
        print L.mt

    ax0.plot(L.mt[0], L.mt[1], 'go')
    showme()


# -----------------------------------
def demo2D2():
    def g(m):
        "rosenbrock"
        d0, d1 = 10. * (m[1] - m[0] ** 2), (1 - m[0])
        return np.array([d0, d1])

    # -----------------------------------
    plt.figure(figsize=(18, 6))
    ax0 = gcf().add_subplot(131)
    ax1 = gcf().add_subplot(132, sharex=ax0, sharey=ax0)


    x = np.linspace(-6., 6., 100)
    y = np.linspace(-6., 6., 110)
    X, Y = np.meshgrid(x, y)
    # Z = np.asarray([g(np.asarray([xx, yy]))[0] for xx, yy in zip(X.flat[:], Y.flat[:])]).reshape(X.shape)
    D0, D1 = l2a(*zip(*[g(np.asarray([xx, yy])) for xx, yy in zip(X.flat[:], Y.flat[:])]))
    D0, D1 = D0.reshape(X.shape), D1.reshape(X.shape)

    ax0.pcolormesh(X, Y, D0, vmin=-abs(D0).max(), vmax=+abs(D0).max(), cmap=plt.cm.spectral)
    gca().set_aspect(1.)
    ax1.pcolormesh(X, Y, D1, vmin=-abs(D1).max(), vmax=+abs(D1).max(), cmap=plt.cm.spectral)
    gca().set_aspect(1.)

    # -----------------------------------
    # data
    d = np.array([0., 0.])  # observed data
    CDinv = np.array([10.0, 1.0]) ** -2.  # observed covariance
    mapr = np.array([4., -4.])  # prior model
    CMinv = np.array([2.0, 2.0]) ** -2.  # coresp. covariance

    # -----------------------------------
    ax0.contour(X, Y, D0, levels=[d[0]], colors="k")
    ax0.contour(X, Y, D0, levels=[d[0] - CDinv[0] ** -.5, d[0] + CDinv[0] ** -.5], colors="k", linestyles="--")
    ax0.contour(X, Y, D1, levels=[d[1]], colors="k", alpha=0.4)
    ax0.contour(X, Y, D1, levels=[d[1] - CDinv[1] ** -.5, d[1] + CDinv[1] ** -.5], colors="k", linestyles="--",
                alpha=0.2)

    ax1.contour(X, Y, D0, levels=[d[0]], colors="k", alpha=0.4)
    ax1.contour(X, Y, D0, levels=[d[0] - CDinv[0] ** -.5, d[0] + CDinv[0] ** -.5], colors="k", linestyles="--",
                alpha=0.2)
    ax1.contour(X, Y, D1, levels=[d[1]], colors="k")
    ax1.contour(X, Y, D1, levels=[d[1] - CDinv[1] ** -.5, d[1] + CDinv[1] ** -.5], colors="k", linestyles="--")

    for ax in ax0, ax1:
        ax.plot(mapr[0], mapr[1], 'ro')
        ax.plot(mapr[0], mapr[1], 'ro')

    ax2 = gcf().add_subplot(133)
    L = LinInv(g=g, d=d, CDinv=CDinv, mapr=mapr, CMinv=CMinv, method=1)
    mt, dt, chi2red = L.upd()
    for n in xrange(30):
        mt, dt, G = L.grad()
        mb, db = L.inv()
        mt, dt, chi2red = L.upd()
        ax0.plot(L.mt[0], L.mt[1], 'w*')
        ax1.plot(L.mt[0], L.mt[1], 'w*')
        ax2.plot(n, L.chi2red, 'ko')
        print L.mt

    ax0.plot(L.mt[0], L.mt[1], 'go')
    ax1.plot(L.mt[0], L.mt[1], 'go')

    ax0.set_title(' $ d_0 = %.2f \pm %.2f $ ' % (d[0], CDinv[0] ** -.5))
    ax1.set_title(' $ d_1 = %.2f \pm %.2f $ ' % (d[1], CDinv[1] ** -.5))
    showme()


# -----------------------------------
def demo1D():
    def g(m):
        """user defined"""
        x = m[0]
        y = 1. - np.exp(-x)
        d = np.asarray([y])
        return d

    # -----------------------------------
    plt.figure(figsize=(12, 6))
    ax0 = gcf().add_subplot(121)

    x = np.linspace(0., 10., 100)
    y = np.asarray([g([xx]) for xx in x])
    gca().plot(x, y)

    # -----------------------------------
    # data
    d = np.array([0.9])
    CDinv = np.array([0.1]) ** -2.
    mapr = np.array([1.e-20])
    CMinv = np.array([1.0]) ** -2.

    ax0.fill_between(minmax(x), (d[0] - CDinv[0] ** -.5) * np.ones(2), (d[0] + CDinv[0] ** -.5) * np.ones(2), color='g',
                     alpha=0.4)
    ax0.fill_betweenx(minmax(y), (mapr[0] - CMinv[0] ** -.5) * np.ones(2), (mapr[0] + CMinv[0] ** -.5) * np.ones(2),
                      color='r', alpha=0.4)
    ax0.plot(minmax(x), d * np.ones(2), 'g')
    ax0.plot(mapr * np.ones(2), minmax(y), 'r')
    ax0.set_xlim(minmax(x))
    ax0.set_ylim(minmax(y))

    # -----------------------------------
    ax1 = gcf().add_subplot(122)
    L = LinInv(g=g, d=d, CDinv=CDinv, mapr=mapr, CMinv=CMinv)
    mt, dt, chi2red = L.upd()
    ax0.plot(L.mt[0], L.dt[0], 'ro')
    for n in xrange(10):
        mt, dt, G = L.grad()
        mb, db = L.inv()
        mt, dt, chi2red = L.upd()
        ax0.plot(L.mt[0], L.dt[0], 'k.')
        ax1.plot(n, L.chi2red, 'ko')

        # -----------------------------------
    ax0.plot(L.mt[0], L.dt[0], 'go')
    showme()


# -----------------------------------
def demogaussreg():
    """show how to fit a gaussian """

    # theory
    class G(object):
        def __init__(self, x):
            self.x = x

        def __call__(self, m):
            a, b, c = m
            return a * np.exp(-0.5 * ((x - b) / abs(c)) ** 2.)

    # -----------------------------------
    truth = [0.463, 0.25, 0.123456]
    x = np.sort(np.random.rand(100))
    g = G(x)
    dat = g(truth)
    nse = 0.05 * np.random.randn(len(dat))
    dat += nse
    sig = 0.05 * np.ones_like(dat)

    # -----------------------------------
    mapr = [0.5, 0.5, 0.5]
    CMinv = np.array([0.01, 0.01, 0.01]) ** -2.
    CDinv = sig ** -2.0

    # -----------------------------------
    gca().plot(x, dat, "ko")
    [gca().plot([xx, xx], [yy - ss, yy + ss], "k_-") for xx, yy, ss in zip(x, dat, sig)]
    # -----------------------------------

    L = LinInv(g=g, d=dat, CDinv=CDinv, mapr=mapr, CMinv=CMinv)
    mt, dt, chi2red = L.upd()
    gca().plot(x, g(L.mt), "r", linewidth=3)
    lastchi2red = np.inf
    while abs(lastchi2red - chi2red) / chi2red > 0.001:  # for n in xrange(100):
        lastchi2red = chi2red
        mt, dt, G = L.grad()
        mb, db = L.inv()
        mt, dt, chi2red = L.upd()
        gca().plot(x, g(L.mt), "k", alpha=0.4)
        print L.mt, chi2red

    gca().plot(x, g(L.mt), "g", linewidth=3)

    showme()


# -----------------------------------
if __name__ == "__main__":
    from tetedenoeud import *


    demogaussreg()
    demo1D()
    demo2D()
    demo2D2()






