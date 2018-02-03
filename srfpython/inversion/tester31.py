from tetedenoeud.utils.display import *
from tetedenoeud.multipro.multipro8 import *
from metropolis2 import *
import numpy as np


"""test Metropolis 2D, parallel"""

####################################################################################
def histogram2d(xflat, yflat, xbins, ybins):
    """quite similar to np.histogram2d with several available normalizations
    """
    if len(xflat) - len(yflat): raise ValueError('')
    H, X, Y = np.histogram2d(x = xflat, y = yflat, bins=(xbins, ybins), normed=True)

    H[np.isnan(H)] = 0.
    return X, Y, H.T

def minmax(X):
    return X.min(), X.max()


x = np.linspace(0.05, 0.3, 200)
y = np.linspace(0.05, 0.4, 215)
X, Y = np.meshgrid(x, y)
#l = LogUniND([0.1, 0.1], [0.3, 0.5], k = 20.)
#l = LogGaussND([0.2, 0.2], [0.025, 1000.], [0.1, 0.1], [0.3, 0.5], k = 20.)
#l = LogGaussNDCov([0.2, 0.2], [0.01, 0.2], [0.1, 0.1], [0.2, 0.5], rho = [[1.0, 0.9], [0.9, 1.0]], k = 20.)
l = LogGaussNDCov(vmeans = [0.2, 0.2],
                  vuncs = [0.1, 0.2],
                  vinfs = [0.1, 0.1],
                  vsups = [0.2, 0.5],
                  rho = [[1.0, 0.9],
                         [0.9, 1.0]],
                  k = 20.,
                  nanbehavior=0)
Z = l.callargs(X, Y)
#ax1 = gcf().add_subplot(131)
#ax1.pcolormesh(X, Y, Z, cmap = plt.cm.jet)#cmap2isocolor(Z, plt.cm.jet))
ax2 = gcf().add_subplot(121)
vmin, vmax = minmax(np.exp(Z))
plt.pcolormesh(X, Y, np.exp(Z), vmin=vmin, vmax=vmax, cmap = plt.cm.jet)
plt.colorbar()


ax3 = gcf().add_subplot(122, sharex = ax2, sharey = ax2)

def gen():
    for nchain in xrange(1):
        M0   = np.random.randn(2)
        MSTD = np.asarray([0.05, 0.05])
        yield Job(nchain, M0, MSTD, nkeep = 15000)

def fun(worker, chainid, M0, MSTD, nkeep):
    def G(model): return np.array([0.])
    def logRHOD(data): return 0.
    logRHOM = l
    models, _, weights, _ = \
        metropolis(M0, MSTD, G, 1, logRHOD, logRHOM,
                   nkeep = nkeep,
                   normallaw = worker.randn,
                   unilaw = worker.rand,
                   chainid = chainid,
                   verbose = True)
    models = models.repeat(weights, axis = 0)

    return models[:, 0], models[:, 1]

X, Y = [], []
with MapAsync(fun, gen(), Nworkers=4, Taskset="0-4") as ma:
    for _, (XX, YY), _, _ in ma:
        X = np.concatenate((X, XX))
        Y = np.concatenate((Y, YY))


X, Y, H = histogram2d(X, Y, x[::4], y[::4])
plt.pcolormesh(X, Y, H, cmap = plt.cm.jet) #vmin=vmin, vmax=vmax,
plt.colorbar()

showme()
