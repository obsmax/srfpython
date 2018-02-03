from tetedenoeud.utils.display import *
from tetedenoeud.multipro.multipro8 import *
from metropolis2 import *
import numpy as np


"""
gaussian regression example : assume incompatible data, and multiple solutions
"""

####################################################################################
def histogram2d(xflat, yflat, xbins, ybins):
    """quite similar to np.histogram2d with several available normalizations
    """
    if len(xflat) - len(yflat): raise ValueError('')
    H, X, Y = np.histogram2d(x = xflat, y = yflat, bins=(xbins, ybins), normed=True)

    H[np.isnan(H)] = 0.
    return X, Y, H.T

#########################################################


#-----------------------
class Theo(object):
    def __init__(self, xpoints):
        self.xpoints = xpoints
    def __call__(self, m):
        #theory
        #d = m[0] * self.xpoints + m[1]
        d = np.exp(-0.5 * ((self.xpoints - m[0]) / m[1]) ** 2.)
        return d
#-----------------------
#data
amin, amax = 0., 1.0
bmin, bmax = 0.05, 0.2
x = np.linspace(0., 1.0, 24)
G = Theo(x)
#-----------------------
y = np.zeros_like(x)
y[:12] = G([0.25, 0.15])[:12] + 0.05 * np.random.randn(len(x[:12]))
y[12:] = G([0.75, 0.08])[12:] + 0.05 * np.random.randn(len(x[12:]))
s = 0.3 * np.ones_like(x)
s[12:] /= 2.
plt.figure(figsize = (12, 6))
ax1 = gcf().add_subplot(121, xlabel = "mean", ylabel = "std"); ax1 = gca()
ax2 = gcf().add_subplot(122, xlabel = "x", ylabel = "y"); ax2 = gca()

for xx,yy,ss in zip(x,y,s):
    ax2.plot([xx, xx], [yy - ss, yy + ss], "_-", color=  [0.5, .5, .5], alpha = 1.0)
ax2.plot(x, y, "wo", markeredgecolor = "k")
showme(0)
#-----------------------
logRHOD = LogGaussND(\
            vmeans = y,
            vuncs  = s,
            vinfs  = y * 0. -10.,
            vsups  = y * 0. +10.)
logRHOM = LogUniND(\
            vinfs  = [amin,  bmin],
            vsups  = [amax,  bmax])
#-----------------------
if False:
    M0 = np.array([np.random.rand() * (amax - amin) + amin, np.random.rand() * (bmax - bmin) + bmin])
    MSTD = np.array([0.1, 0.1])
    nkeep = 10000
    models, datas, weights, llks = metropolis(M0, MSTD, G, len(G.xpoints), logRHOD, logRHOM, nkeep = nkeep)
    models = models.repeat(weights, axis = 0)
    datas  = datas.repeat(weights, axis = 0)
    llks   = llks.repeat(weights)
    A = models[:, 0]
    B = models[:, 1]
else: #parallel
    def gen():
        for nchain in xrange(48):
            M0 = np.array([np.random.rand() * (amax - amin) + amin, np.random.rand() * (bmax - bmin) + bmin])
            MSTD = np.array([0.1, 0.1])
            nkeep = 1000#00
            yield Job(nchain, M0, MSTD, nkeep = nkeep)

    def fun(worker, chainid, M0, MSTD, nkeep):
        G = Theo(x)
        logRHOD = LogGaussND(\
                    vmeans = y,
                    vuncs  = s,
                    vinfs  = y * 0. -10.,
                    vsups  = y * 0. +10.)
        logRHOM = LogUniND(\
                    vinfs  = [amin,  bmin],
                    vsups  = [amax,  bmax])
        models, datas, weights, llks = metropolis(M0, MSTD, G, len(G.xpoints), logRHOD, logRHOM, nkeep = nkeep, normallaw = worker.randn, unilaw = worker.rand, chainid = chainid)
        models = models.repeat(weights, axis = 0)
        datas  = datas.repeat(weights, axis = 0)
        llks   = llks.repeat(weights)
        return models, datas, weights, llks

    models = None
    datas  = None
    with MapAsync(fun, gen(), Nworkers=4) as ma:
        for _, (ms, ds, ws, ls), _, _ in ma:
            if models is None:
                models, datas = ms, ds
            else:
                models = np.concatenate((models, ms), axis = 0)
                datas = np.concatenate((datas, ds), axis = 0)
    A = models[:, 0]
    B = models[:, 1]
#----------------------- data histogram
xbins = np.linspace(x.min(), x.max(), 100)
ybins = np.linspace((y - s).min(), (y + s).max(), 110)
Xbins, Ybins = np.meshgrid(xbins, ybins)
Z = np.zeros_like(Xbins)
J = np.arange(len(xbins))
G1 = Theo(xbins)
for a, b in zip(A, B):
    y = G1([a, b])
    I = np.argmin(abs(y - Ybins), axis = 0)
    Z[I, J] += 1
cmap = plt.cm.CMRmap #cmap2isocolor(Z, cmap = plt.cm.CMRmap)
ax2.pcolormesh(xbins, ybins, Z, cmap = cmap)

#----------------------- model histogram
X, Y, H = histogram2d(A, B, np.linspace(amin, amax, 100), np.linspace(bmin, bmax, 100))
#cmap = cmap2isocolor(H, cmap = plt.cm.CMRmap)
cmap = plt.cm.CMRmap
ax1.pcolormesh(X, Y, H, cmap = cmap)

showme()
