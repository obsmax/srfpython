from tetedenoeud.display import *
from metropolis2 import *


"""
test Metropolis algorithm in 1D
"""

#l = LogGauss(vmean = 0., vunc = 1.0, vinf = -0.1, vsup = 3.0, k = 1000.0, nanbehavior=  0) #not adapted for metropolis
if 1:
    l = LogGaussND(vmeans=[0.], vuncs=[0.25],
               vinfs=[-0.1], vsups=[3.0],
               k=10000.0, nanbehavior=0)
    logRHOM = l
elif 1:
    l = LogUniND(vinfs=[-0.1], vsups=[1.],
                 k=10000.0, nanbehavior=0)
    logRHOM = l


x = np.linspace(-.5, 3., 1000)
y = np.exp([logRHOM([xx]) for xx in x])

A = np.sum(y * (x[1] - x[0]))
gca().plot(x, y / A)


def G(model):
    return np.array([0.])

def logRHOD(data):
    return 0.


M0 = np.array([0.])
MSTD = np.array([20.])
start = time.time()
models, _, weights, _ = metropolis(M0, MSTD, G, 1,
                                   logRHOD, logRHOM,
                                   nkeep = 500000,
                                   HL = 300, IK0 = 0.65,
                                   adjustspeed=0.1,
                                   MPMIN = .1, MPMAX = 10.0)
print "total time %f" % (time.time() -start)
models = models.repeat(weights, axis = 0)
hist, bin_edges = np.histogram(models[:, 0], bins = 30, density=True)
for mid, width, height in zip(.5 * (bin_edges[:-1] + bin_edges[1:]), bin_edges[1:]-bin_edges[:-1], hist):
    gca().bar(mid, height, width, bottom=0, align='center', color = "k", alpha = .4)



gca().set_xlim(-.5, 2.)
showme()
