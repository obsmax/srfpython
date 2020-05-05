import numpy as np
from scipy.sparse import diags, csr_matrix, csc_matrix, \
    save_npz as save_sparse_npz, load_npz as load_sparse_npz
from srfpython.standalone.multipro8 import Job, MapSync
from srfpython.HerrMet.paramfile import load_paramfile
from srfpython.HerrMet.datacoders import makedatacoder, Datacoder_log
from srfpython.HerrMet.theory import Theory

"""
compute sensitivity kernels
"""


paramterizer_strings = np.load('paramterizer_strings.npy')
datacoder_strings = np.load('datacoder_strings.npy')
parameterizers = [load_paramfile(paramterizer_string, verbose=False)[0]
                  for paramterizer_string in paramterizer_strings]
datacoders = [makedatacoder(datacoder_string, Datacoder_log)
              for datacoder_string in datacoder_strings]
theorys = [Theory(parameterizer=parameterizer, datacoder=datacoder)
           for parameterizer, datacoder in zip(parameterizers, datacoders)]


nds = np.asarray([len(dc.values) for dc in datacoders], int)
nps = np.asarray([len(p.MMEAN) for p in parameterizers], int)
G_shape = (nds.sum(), nps.sum())
G_rows = np.hstack((0, nds.cumsum()[:-1]))  # index of the first row of G for each node
G_cols = np.hstack((0, nps.cumsum()[:-1]))  # index of the first col of G for each node

# =============================
def gen():
    for nnode, (parameterizer, datacoder, theory) in enumerate(zip(parameterizers, datacoders, theorys)):
        m = parameterizer.MMEAN
        # if True:
        #     gm = theory(m)
        #     print(gm)
        #     if np.isnan(gm).any():
        #         dm = parameterizer.inv_to_depthmodel(m)
        #         laws = datacoder.inv_to_laws(gm)
        #         import matplotlib.pyplot as plt
        #         plt.figure()
        #         dm.show(plt.gca())
        #         plt.figure()
        #         for law in laws:
        #             law.show(plt.gca(), period=True)
        #         plt.show()


        yield Job(nnode, theory, m, G_rows[nnode], G_cols[nnode])


def fun(nnode, theory, m, G_first_row, G_first_col):
    # compute Frechet derivatives near parameterizer.MMEAN
    fd = theory.frechet_derivatives(m=m, gm=None)
    fd_col_ind, fd_row_ind = np.meshgrid(range(fd.shape[1]), range(fd.shape[0]))

    # store indices and values to fill the sparse matrix G
    G_rows = G_first_row + fd_row_ind.flat[:]
    G_cols = G_first_col + fd_col_ind.flat[:]
    G_datas = fd.flat[:]
    return nnode, G_rows, G_cols, G_datas

# [fun(*_.args, **_.kwargs) for _ in list(gen())]

G_row_ind = np.array([], int)
G_col_ind = np.array([], int)
G_fd_data = np.array([], float)
with MapSync(fun, gen()) as ma:
    for jobid, (nnode, G_rows, G_cols, G_datas), _, _ in ma:
        G_row_ind = np.concatenate((G_row_ind, G_rows))
        G_col_ind = np.concatenate((G_col_ind, G_cols))
        G_fd_data = np.concatenate((G_fd_data, G_datas))


G = csr_matrix((G_fd_data, (G_row_ind, G_col_ind)), shape=G_shape)
save_sparse_npz('G.npz', G)

if False:
    assert raw_input('sure?') == "y"
    import matplotlib.pyplot as plt
    plt.figure()
    #plt.colorbar(plt.imshow(G.toarray(), vmin=-1.5, vmax=1.5, cmap=plt.get_cmap('RdBu')))
    # plt.plot(G.sum(axis=1))
    G = G.toarray()
    plt.colorbar(plt.imshow(G))
    plt.figure()
    plt.plot(G.sum(axis=1))
    plt.plot(G.sum(axis=0))
    plt.show()
    exit()

