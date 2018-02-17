from srfpython.standalone.display import *
import numpy as np
from metropolis2 import *

"""
test the 1D LogUni and LogGauss pdfs 
"""

v = np.linspace(-1., 1., 10000)

if 0:
    #nans expected
    l0 = LogUni(vinf=0., vsup=0.1, k=10., nanbehavior=2)
    l1 = LogUni(vinf=0.1, vsup=0.2, k=30., nanbehavior=2)
    l2 = LogGauss(vmean=0.22, vunc=0.12, vinf=0.2, vsup=0.3, k=30., nanbehavior=2)

    V = np.copy(v)
    V[-1] = np.nan
else:
    #no nans expected
    l0 = LogUni(vinf=0., vsup=0.1, k=10., nanbehavior= 0)
    l1 = LogUni(vinf=0.1, vsup=0.2, k=30., nanbehavior= 0)
    l2 = LogGauss(vmean=0.22, vunc=0.12, vinf=0.2, vsup=0.3, k=30., nanbehavior= 0)
    V = np.copy(v)


print "l0", l0
print "l1", l1
print "l2", l2

gcf().add_subplot(211, ylabel="log(pdf)")
gca().plot(v, l0.calln(V), label="l0 (uni, k=10)")
gca().plot(v, l1.calln(V), label="l1 (uni, k=30)")
gca().plot(v, l2.calln(V), label="l2 (trunc. gauss., k=30)")
gca().set_ylim(-100., 2.)
gca().legend()
gcf().add_subplot(212, sharex = gca(), ylabel="pdf")
gca().plot(v, np.exp(l0.calln(V)), label="l0 (uni, k=10)")
gca().plot(v, np.exp(l1.calln(V)), label="l1 (uni, k=30)")
gca().plot(v, np.exp(l2.calln(V)), label="l2 (trunc. gauss., k=30)")
gca().legend()
gca().set_xlim(-.05, .4)
showme()
