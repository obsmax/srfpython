import numpy as np
from scipy import sparse

"""
test for matrix products
"""


A = np.array([[1, 2, 3], [0, 4, 5], [0, 0, 6]])

j, i = [_.flat[:]
        for _ in np.meshgrid(np.arange(A.shape[0]),
                             np.arange(A.shape[1]))]

I = A.flat[:] != 0
A = sparse.csc_matrix((A.flat[I], (i[I], j[I])), shape=(3, 3))
print(sparse.issparse(A))
x = np.array([1, 2, 3])

# xTA = np.dot(x, A) # NO !!!!
# xTA = np.matmul(x, A) # no
# xTA = np.dot(x, A.toarray())  # yes
# xTA = x * A  # yes
# xTA = A.__rmul__(x)  # yes
# xTA = x.__mul__(A)  # no
xTA = A.__rmatmul__(x)  # yes

print(type(xTA), xTA.__class__.__name__, xTA)
