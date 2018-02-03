from tetedenoeud.utils.display import *
from tetedenoeud.multipro.multipro8 import *
from metropolis2 import *
import numpy as np


"""
test for failing theory or nans in data array
use the Programming Error exception to test the behavior of the code
use failure test or happy failure to mimic unsuited models
"""
nofail = False #if True, the chain will continue running even though programming error is suspected
nanbehavior = 2 #recommended : means that models or datas with nans are penalized according to the number of nans
                #0 means that nans are producing errors, so the metropolis algo cannot compare models since they are all equally bad


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
        if False:
            raise Exception('Programming Error test')

        if m[1] >= 0.16:
            return np.nan * np.zeros_like(self.xpoints)
            raise Exception('Theory failure test')
        elif m[1] <= 0.08:
            return np.nan * np.zeros_like(self.xpoints)
        elif np.random.rand() < 0.01:
            return np.nan * np.zeros_like(self.xpoints) #raise Exception('Happy failure test : the failure that occurs just for fun')
        d = np.exp(-0.5 * ((self.xpoints - m[0]) / m[1]) ** 2.)
        return d
#-----------------------
#data
amin, amax = 0., 1.0
bmin, bmax = 0.05, 0.2
x = np.linspace(0., 1.0, 12)
G = Theo(x)
#-----------------------
y = np.zeros_like(x)
y = np.exp(-0.5 * ((x - 0.5) / 0.1) ** 2.) + 0.1 * np.random.randn(len(x))
s = 0.3 * np.ones_like(x)
plt.figure(figsize = (12, 6))
ax1 = gcf().add_subplot(121, xlabel = "mean", ylabel = "std"); ax1 = gca()
ax2 = gcf().add_subplot(122, xlabel = "x", ylabel = "y"); ax2 = gca()

for xx,yy,ss in zip(x,y,s):
    ax2.plot([xx, xx], [yy - ss, yy + ss], "_-", color=  [0.5, .5, .5], alpha = 1.0)
ax2.plot(x, y, "wo", markeredgecolor = "k")
showme(0)
#-----------------------
if True:
    def gen():
        for nchain in xrange(12):
            M0 = np.array([np.random.rand() * (amax - amin) + amin, np.random.rand() * (bmax - bmin) + bmin])
            MSTD = np.array([0.1, 0.1])
            nkeep = 1000
            yield Job(nchain, M0, MSTD, nkeep = nkeep)

    def fun(worker, chainid, M0, MSTD, nkeep):
        G = Theo(x)
        logRHOD = LogGaussND(\
                    vmeans = y,
                    vuncs  = s,
                    vinfs  = y * 0. -10.,
                    vsups  = y * 0. +10.,
                    nanbehavior = nanbehavior)
        logRHOM = LogUniND(\
                    vinfs  = [amin,  bmin],
                    vsups  = [amax,  bmax],
                    nanbehavior = nanbehavior)
        models, datas, weights, llks = metropolis(M0, MSTD, G, len(G.xpoints), logRHOD, logRHOM, nkeep = nkeep, normallaw = worker.randn, unilaw = worker.rand, chainid = chainid, nofail = nofail)
        models = models.repeat(weights, axis = 0)
        datas  = datas.repeat(weights, axis = 0)
        llks   = llks.repeat(weights)
        return models, datas, weights, llks

    models = None
    datas  = None
    with MapAsync(fun, gen()) as ma:
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
#G1 = Theo(xbins)
for a, b in zip(A, B):
    try: y = np.exp(-0.5 * ((xbins - a) / b) ** 2.)#G1([a, b])#G1([a, b])
    except KeyboardInterrupt: raise
    except: continue
    I = np.argmin(abs(y - Ybins), axis = 0)
    Z[I, J] += 1

cmap = plt.cm.CMRmap #cmap2isocolor(Z, cmap = plt.cm.CMRmap)
ax2.pcolormesh(xbins, ybins, Z, cmap = cmap)
#----------------------- model histogram
X, Y, H = histogram2d(A, B, np.linspace(amin - 0.33, amax + 0.33, 100), np.linspace(bmin - 0.08, bmax + 0.08, 100))#, normalization="pdf")
cmap = plt.cm.CMRmap #cmap2isocolor(H, cmap = plt.cm.CMRmap)
ax1.pcolormesh(X, Y, H, cmap = cmap)
showme()

