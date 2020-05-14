import numpy as np
from scipy import sparse as sp
from scipy.sparse import linalg as splinalg
from srfpython.standalone.multipro8 import *

# not omptimal : need to recmpute the same eigen values many times


def concurent_eig_search(A):

    def search(sign):
        # A is readonly in this function
        valps = []
        k = 1
        while k < A.shape[1]:
            if sign > 0:
                valp = splinalg.eigsh(A, which='LA', k=k, return_eigenvectors=False)[0]
            else:
                valp = splinalg.eigsh(A, which='SA', k=k, return_eigenvectors=False)[0]

            if valp * sign > 0:
                valps.append(valp)
            else:
                break
            k += 1
        return sign, valps

    def job_generator():
        yield Job(1)   # one thread searches for positive eigenvalues
        yield Job(-1)  # the other searches for negative eigenvalues

    with MapAsync(search, job_generator(), Nworkers=2) as ma:
        for jobid, (sign, valps), _, _ in ma:
            break  # the thread that finishes first has won the battle

    return sign, valps


if __name__ == '__main__':

    # CM = sp.diags(np.arange(-6, 10., 1))
    CM = sp.csc_matrix(np.random.randn(10, 10))
    print(concurent_eig_search(CM))
