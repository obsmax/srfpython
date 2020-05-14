import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse as sp
from scipy.sparse import linalg as splinalg

a = 1.0
b = 0.25
c = 0.01

d = 0.8
e = 0.5
f = 0.01

o = 0.01

C = np.array([
    [a,   b,   c,        d,   e,   f,        o,   o,   o],
    [b,   a,   b,        e,   d,   e,        o,   o,   o],
    [c,   b,   a,        f,   e,   d,        o,   o,   o],

    [d,   e,   f,        a,   b,   c,        d,   e,   f],
    [e,   d,   e,        b,   a,   b,        e,   d,   e],
    [f,   e,   d,        c,   b,   a,        f,   e,   d],

    [o,   o,   o,        d,   e,   f,        a,   b,   c],
    [o,   o,   o,        e,   d,   e,        b,   a,   b],
    [o,   o,   o,        f,   e,   d,        c,   b,   a]
    ])


valp, vectp = np.linalg.eigh(C)
I = np.argsort(valp)[::-1]
valp = valp[I]
vectp = vectp[:, I]

plt.plot(valp)
plt.show()
print(valp)
print(vectp)

Ipos = valp > 0
Cfix = np.asarray(
    np.asmatrix(vectp[:, Ipos]) * \
    np.asmatrix(np.diag(valp[Ipos])) * \
    np.asmatrix(vectp[:, Ipos]).T)
plt.figure()
plt.subplot(121)
plt.colorbar(plt.imshow(C))
plt.subplot(122, sharex=plt.gca(), sharey=plt.gca())
plt.colorbar(plt.imshow(Cfix))
plt.show()


exit()













C = np.array([
    [1.,   1.,   0.],
    [1.,   1.,   1.],
    [0.,   1.,   1.]])



C = sp.csc_matrix(C)
C.prune()
valp, vectp = splinalg.eigsh(C, k=1, which='SA')
print(valp, vectp)


exit()
valp, vectp = np.linalg.eigh(C)
I = np.argsort(valp)[::-1]
valp = valp[I]
vectp = vectp[:, I]

plt.plot(valp)
plt.show()
exit()



r2 = np.sqrt(2.) / 2.

V = np.array([[r2, -r2],
              [r2,  r2]])
V = np.asmatrix(V)

L = np.array([[1., 0.],
              [0., 0.1]])

L = np.asmatrix(L)

C = V * L * V.T
print(C)


Cinv = np.linalg.inv(C)
print(Cinv)


x = np.linspace(-5., 5., 100)
y = np.linspace(-5., 5., 110)

x, y = np.meshgrid(x, y)

X = np.column_stack((x.flat[:], y.flat[:]))
X = np.asmatrix(X)
print(X.shape)


K = (Cinv * X.T).T

K = np.asarray(K)
X = np.asarray(X)
f = np.exp(-.5 * ((K * X).sum(axis=1)))
f = f.reshape(x.shape)

plt.gca().set_aspect(1.0)
plt.contourf(x, y, f)
plt.show()
