from tetedenoeud.utils.display import *
from tetedenoeud.multipro.multipro8 import *
from metropolis2 import *
import numpy as np


"""add constraint on secondary parameters, 2D
   model = x, y
   target = (x ** 2. + y ** 2.) ** 0.5 = 1.0 +- 0.1
   constraint = abs(x - y) follows a uniform law between two values

   > a model is [x, y]
   > the data is G(model) is [d]
   > the prior pdf on the model is bypassed so that additional constraints are added to some relations
     between the model parameters

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


# data pdf
def G(model):
    x, y = model
    d = np.sqrt(x ** 2. + y ** 2.)
    return np.array([d])


logRHOD = LogGaussND(vmeans=[1.0],
                     vuncs=[0.1],
                     vinfs=[0.5],
                     vsups=[1.5],
                     k=1000.0,
                     nanbehavior=0)

# model pdf
# prior constraint on x      y    abs(x-y)
l = LogUniND(vinfs=[-10., -10.0, 0.5],
             vsups=[10., 10.0, 1.5],
             k=1000.0,
             nanbehavior=0)

def logRHOM(model):
    x, y = model
    return l([x, y, np.abs(x - y)])


if True:
    def gen():
        for nchain in xrange(4):
            M0   = np.random.randn(2)
            MSTD = np.asarray([1., 1.])
            yield Job(nchain, M0, MSTD)

    def fun(worker, chainid, M0, MSTD):
        models, _, weights, _ = \
            metropolis(M0, MSTD, G, 1, logRHOD, logRHOM,
                       nkeep=10000,
                       HL=1000,
                       IK0=0.25,
                       MPMIN=0.001,
                       MPMAX=1000.0,
                       adjustspeed=0.1,
                       debug=False,
                       normallaw = worker.randn,
                       unilaw = worker.rand,
                       chainid = chainid,
                       verbose = True)
        models = models.repeat(weights, axis = 0)

        return models

models=None
with MapAsync(fun, gen()) as ma:
    for _, m, _, _ in ma:
        if models is None:
            models = m
        else:
            models = np.concatenate((models, m), axis=0)

X, Y, H = histogram2d(models[:, 0], models[:, 1], 100, 100)#, normalization="pdf")
cmap = plt.cm.CMRmap #cmap2isocolor(H, cmap = plt.cm.CMRmap)
gca().pcolormesh(X, Y, H, cmap = cmap)
gca().set_aspect(1.0)
showme()

