from tetedenoeud.display import *
import numpy as np
from metropolis2 import *

"""
test the 1D LogUni and LogGauss pdfs 
"""

v = np.linspace(-1., 1., 10000)

if 0:
    #nans expected
    l0 = LogUni(0., 0.1, k = 20., nanbehavior= 2)
    l1 = LogGauss(0.05, 100000., 0., 0.1, k = 20., nanbehavior= 2)
    l2 = LogGauss(0., 0.1, 0., 0.1, k = 20., nanbehavior= 2)

    V = np.copy(v)
    V[-1] = np.nan
else:
    #no nans expected
    k = 30.
    nanbehavior = 0
    l0 = LogUni(0., 0.1, k = k, nanbehavior= 0)
    l1 = LogGauss(0.05, 100000., 0., 0.1, k = k, nanbehavior= 0)
    l2 = LogGauss(0., 0.012, 0., 0.1, k = k, nanbehavior= 0)
    V = np.copy(v)

gcf().add_subplot(211)
gca().plot(v, l0.calln(V), label = "l0")
gca().plot(v, l1.calln(V), label = "l1")
gca().plot(v, l2.calln(V), label = "l2")
gca().set_ylim(-5., 2.)
gca().legend()
gcf().add_subplot(212, sharex = gca())
gca().plot(v, np.exp(l0.calln(V)), label = "l0")
gca().plot(v, np.exp(l1.calln(V)), label = "l1")
gca().plot(v, np.exp(l2.calln(V)), label = "l2")
gca().legend()
gca().set_xlim(-.02, .12)
showme()
