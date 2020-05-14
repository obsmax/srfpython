import numpy as np
from srfpython.standalone.multipro8 import *

# shared variable (read only)
bigone = np.ones(1024 ** 3 // 8, float)   # 1Gb array


def jobgen():
    for i in range(10):
        yield Job(i)


def fun(i):
    # use bigone in readonly without passing it to fun as input
    # => the array is shared between all processes

    # this does not work if bigone is modified inside this function

    print(">", bigone[:10])
    start = time.time()
    while time.time() - start < 2.:
        0. + 0.
    return i


with MapAsync(fun, jobgen(), Nworkers=10) as ma:
    list(ma)


