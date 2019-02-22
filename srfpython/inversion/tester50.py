from srfpython import *
from neldermead2 import *
import numpy as np


def gauss2D2(x, y, xo, yo, sigmaj, sigmin, alphadeg):
    """plus pratique, pas de coeff de correl
    xo, yo = max of the gaussian
    sigmaj = standard deviation in the major eigenvector direction
    sigmin = standard deviation in the minor eigenvector direction
    alphadeg = angle in degrees, counterclockwise from the positive x axis
    """
    assert sigmaj >= sigmin

    deg2rad = np.pi / 180.
    x = x - xo
    y = y - yo
    c =  np.cos(alphadeg * deg2rad)
    s =  np.sin(alphadeg * deg2rad)
    xp = (y * c - x * s) / sigmin
    yp = (x * c + y * s) / sigmaj
    return np.exp(-0.5 * (xp ** 2. + yp ** 2.))


class fun2D():
    def __init__(self, N=20):
        if False:
            self.N = N
            self.As = np.random.randn(self.N)
            self.xos, self.yos = zip(*[np.random.rand(2) * 20. - 10. for _ in xrange(self.N)])
            self.sigmamins, self.sigmajs = zip(*[np.sort(np.random.rand(2) * 3. + .1) for _ in xrange(self.N)])
            self.alphadegs = np.random.rand(self.N) * 360.
        else:
            self.As = np.array([4., 1.])
            self.xos = np.array([+8.5, +5.0])
            self.yos = np.array([-2.5, -5.0])
            self.sigmamins = np.array([.1, 4.])
            self.sigmajs   = np.array([3., 4.])
            self.alphadegs = np.array([35., 0.])
            self.N =  len(self.As)

    def __call__(self, m):
        x, y = m
        z = np.zeros_like(x)
        for n in xrange(self.N):
            z += self.As[n] * gauss2D2(x, y, self.xos[n], self.yos[n], self.sigmajs[n], self.sigmamins[n],
                                       self.alphadegs[n])

        z += 0. * np.cos(x*10.) * np.sin(y*10.)
        return z


class RosenBrock(object):
    "rosenbrock * -1"
    def __init__(self, a=1, b=100):
        self.a, self.b = a, b

    def __call__(self, m):
        x, y = m
        return -1 * ((self.a - x) ** 2. + self.b * (y - x ** 2.) ** 2.)


def test2D():
    x = np.linspace(-2.5, 10., 300)
    y = np.linspace(-10., 2.5, 300)
    X, Y = np.meshgrid(x, y)
    if True:
        f = fun2D()
    else:
        f = RosenBrock()
    Z = f(np.asarray([X, Y]))

    plt.figure(figsize=(12, 6))
    ax0 = gcf().add_subplot(121);
    ax0 = gca()
    ax1 = gcf().add_subplot(122)

    ax0.pcolormesh(X, Y, Z, vmin=Z.min(), vmax=Z.max(), cmap=plt.cm.spectral)
    #ax0.contour(X, Y, Z, vmin=Z.min(), vmax=Z.max(), colors="gray")

    ax0.set_aspect(1.0)

    def G(m):
        return 0. # fake theory
    def logrhod(d):
        return 0. # no data misfit
    def logrhom(m):
        return f(m) # fake prior pdf = the function to maximize

    plt.sca(ax0)
    models, datas, llks = neldermead(np.asarray([0., 0.]),
                   DM=1.,
                   G=G,
                   ND=1,
                   logRHOD=logrhod,
                   logRHOM=logrhom,
                   alpha=1.5,
                   beta=0.5,
                   gamma=2.0,
                   interrupt=1e-3)

    ax0.plot(models[:, 0], models[:, 1], 'w')
    ax1.plot(llks, "ks")
    """
    print models
    print datas
    print llks

    for n, (Mi, llki, (Mis, Dis, llks)) in enumerate(g):
        xtri = np.concatenate((Mis[0, :], [Mis[0, 0]]))
        ytri = np.concatenate((Mis[1, :], [Mis[1, 0]]))
        ax0.plot(xtri, ytri, "w")
        ax0.plot(Mi[0], Mi[1], "k.")
        ax1.plot(n, llki, "ko")

    print Mi[0], Mi[1]
    ax0.plot(Mi[0], Mi[1], "g*")
    ax1.plot(n, llki, "g*")
    ax1.grid(True)
    """

    showme()


test2D()
