import numpy as np
import matplotlib.pyplot as plt

"""
comment on Tarantola 2005, section 3.4.5 p80
on the inversion of a large matrix
"""


def Ainv_dot_b(A, b):
    """
    computes A^-1 * b without inverting A
    for numerical stability
    see Tarantola 2005, section 3.4.5 p80
    :param A:
    :param b:
    :return:
    """
    return np.linalg.solve(a=A, b=b)

# ========== construct a non diagonal covariance matrix
n = 1000
z = np.linspace(0., 3., n)
sigma = np.linspace(0.1, 0.1, n)
l = 1.0

sigmaT = sigma[:, np.newaxis]
zT = z[:, np.newaxis]

if l > 0:
    #rho = np.exp(-0.5 * ((z - zT) / l) ** 2.0)
    rho = np.exp(-np.abs(z - zT) / l)
elif l == 0:
    rho = np.diag(np.ones(n))
else:
    raise ValueError

C = sigma * sigmaT * rho

det = np.linalg.det(C)
# ========== compute the determinant
print(det, det == 0)  # => 0, True
Cinv = np.linalg.inv(C)  # => inaccurate
m = np.ones(n)
print (np.dot(m, np.dot(Cinv, m)))  # => negative ???
print (np.dot(m, Ainv_dot_b(C, m)))   # => looks ok

