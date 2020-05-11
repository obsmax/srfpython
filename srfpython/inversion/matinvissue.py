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


# def Ainv_dot_B(A, B):
#     ans = np.zeros((A.shape[0], B.shape[1]))
#     for j in range(B.shape[1]):
#         ans[:, j] = Ainv_dot_b(A, b=B[:, j])
#     return ans


# ========== construct a non diagonal covariance matrix
n = 1000
z = np.linspace(0., 3., n)
sigma = np.linspace(0.1, 0.1, n)
l = 1.0

sigmaT = sigma[:, np.newaxis]
zT = z[:, np.newaxis]

if l > 0:
    # rho = np.exp(-0.5 * ((z - zT) / l) ** 2.0)
    rho = np.exp(-np.abs(z - zT) / l)
elif l == 0:
    rho = np.diag(np.ones(n))
else:
    raise ValueError

C = sigma * sigmaT * rho
# L = np.linalg.cholesky(C)
# plt.figure()
# plt.subplot(121)
# plt.imshow(C)
# plt.subplot(122)
# plt.imshow(L)
# plt.show()
#
#
# det = np.linalg.det(C)
# ========== compute the determinant
# print(det, det == 0)  # => 0, True
Cinv = np.linalg.inv(C)  # => inaccurate
m = np.random.rand(n)
# m = np.ones(n)

print (np.dot(m, np.dot(Cinv, m)))  # => negative ???
print (np.dot(m, Ainv_dot_b(C, m)))  # => looks ok

x = Ainv_dot_b(C, m)
plt.figure()
plt.subplot(211)
plt.plot(m)
plt.plot(np.dot(C, x))
plt.subplot(212, sharex=plt.gca())
plt.plot(m - np.dot(C, x))
plt.show()