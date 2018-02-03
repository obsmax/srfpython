from tetedenoeud.utils.display import *
from tetedenoeud.multipro.multipro8 import *
from metropolis2 import *
import numpy as np

"""add constraint on secondary parameters
x is the parameter to invert
we impose a constraint on x and one on 2*x
"""

if False:
    # prior constraint on x only
    l = LogUni(vinf=0.,
               vsup=2.5,
               k=1000.0,
               nanbehavior=0)


    def logRHOM(model):
        x = model[0]
        return l(x)

else:
    #prior constraint on x    2*x
    l = LogUniND(vinfs=[-10., 0.0],
                 vsups=[10.0, 5.0], #=> x should be only between 0. and 2.5
                 k=1000.0,
                 nanbehavior=0)


    def logRHOM(model):
        x = model[0]
        y = 2. * x
        return l([x, y])


def G(model):
    return np.array([0.])


def logRHOD(data):
    return 0.


x = np.linspace(-10., 10., 1000)
y = np.exp([logRHOM([xx]) for xx in x])

A = np.sum(y * (x[1] - x[0]))
gca().plot(x, y / A)


M0 = np.array([0.])
MSTD = np.array([10.])
models, _, weights, _ = \
    metropolis(M0, MSTD, G, 1, logRHOD, logRHOM,
               nkeep=50000,
               HL=100,
               IK0=0.65,
               MPMIN=0.001,
               MPMAX=1000.0,
               adjustspeed=0.5)
models = models.repeat(weights, axis=0)
#p = PDF(models[:, 0], 200)
#p.plot(gca(), alpha = 0.5)

hist, bin_edges = np.histogram(models[:, 0], bins =30, density=True)
for mid, width, height in zip(.5 * (bin_edges[:-1] + bin_edges[1:]), bin_edges[1:]-bin_edges[:-1], hist):
    gca().bar(mid, height, width, bottom=0, align='center', color = "k", alpha = .4)


# compare with the distribution that would be obtained with random.rand
N = models.shape[0]
models = np.random.rand(N).reshape((N, 1)) * 2.5 + 3.
hist, bin_edges = np.histogram(models[:, 0], bins =30, density=True)
for mid, width, height in zip(.5 * (bin_edges[:-1] + bin_edges[1:]), bin_edges[1:]-bin_edges[:-1], hist):
    gca().bar(mid, height, width, bottom=0, align='center', color = "g", alpha = .4)


showme()
