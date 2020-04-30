from __future__ import print_function
from builtins import  input
import sys, glob, os
import numpy as np
from srfpython.HerrMet.nodefile import NodeFileString, NodeFile
from srfpython.HerrMet.theory import Theory
from srfpython.HerrMet.paramfile import load_paramfile
from srfpython.HerrMet.datacoders import makedatacoder, Datacoder_log
from srfpython.standalone.asciifile import AsciiFile
from srfpython.HerrMet.files import HERRMETEXTRACTPDFMODELFILE, HERRMETTARGETFILE, ROOTNAME
from srfpython.depthdisp.depthmodels import depthmodel_from_mod96, brocher2005, depthmodel_from_arrays
from scipy.sparse import csr_matrix, save_npz as save_sparse_npz, load_npz as load_sparse_npz

"""
"""
rootdir = "/home/max/prog/git/srfpython/tutorials/02_cube_inversion_example/inversion"

# search the extracted files
search_path = \
    HERRMETEXTRACTPDFMODELFILE.format(
    rootname=os.path.join(rootdir, ROOTNAME.format(node="*")),
    extract_mode="best",
    extract_limit=1000,
    extract_llkmin=0,
    extract_step=1,
    percentile=0.5)

ztop = np.linspace(0., 3.0, 30)

extract_files = glob.glob(search_path) #[:5]
current_row = 0
current_col = 0
row_ind = np.array([], int)
col_ind = np.array([], int)
fd_data = np.array([], float)
for nnode, m96 in enumerate(extract_files):

    # find the corresponding target file
    s96 = HERRMETTARGETFILE.format(
        rootname=os.path.dirname(m96))
    s96 = os.path.join(rootdir, s96)
    assert os.path.isfile(s96)

    # load the extracted file
    dm = depthmodel_from_mod96(m96)

    # resample to ztop
    vs = dm.vs.interp(ztop)

    # # force vp=f(vs), rh=f(vp(vs)) for now
    # vp, rh = brocher2005(vs)
    #
    # # recreate a depthmodel
    # dm = depthmodel_from_arrays(ztop, vp, vs, rh)

    # write a parameter file for the optimizer
    parameter_string = """
    #met NLAYER = {}
    #met TYPE = 'mZVSVPvsRHvp'
    #met VPvs = 'lambda VS: 0.9409 + 2.0947 * VS - 0.8206 * VS ** 2.0 + 0.2683 * VS ** 3.0 - 0.0251 * VS ** 4.0'
    #met RHvp = 'lambda VP: 1.6612 * VP - 0.4721 * VP ** 2.0 + 0.0671 * VP ** 3.0 - 0.0043 * VP ** 4.0 + 0.000106 * VP ** 5.0'
    #fld KEY     VINF          VSUP
    #unt []      []            []
    #fmt %s      %f            %f
    """.format(len(ztop)).replace('    #', '#')

    for i in range(1, len(ztop)):
        parameter_string += "-Z{} {} {}\n".format(i, -ztop[i], -ztop[i])  # add locked depth interfaces

    for i in range(len(ztop)):
        parameter_string += "VS{} {} {}\n".format(i, vs[i]-0.01, vs[i]+0.01)  # so that MMEAN corresponds to the extracted vs

    # initiate the parameterizer and datacoder
    parameterizer, _ = load_paramfile(parameter_string, verbose=False)
    datacoder = makedatacoder(s96, which=Datacoder_log)

    # initiate the theory operator (g)
    theory = Theory(parameterizer=parameterizer, datacoder=datacoder)

    # compute Frechet derivatives
    fd = theory.frechet_derivatives(m=parameterizer.MMEAN, gm=None)
    fd_col_ind, fd_row_ind = np.meshgrid(range(fd.shape[1]), range(fd.shape[0]))

    # qc
    if 0:
        import matplotlib.pyplot as plt
        plt.gcf().clf()
        plt.ion()
        plt.imshow(fd)
        plt.show()
        input('pause')

    # fill the sparse matrix
    row_ind = np.concatenate((row_ind, current_row + fd_row_ind.flat[:]))
    col_ind = np.concatenate((col_ind, current_col + fd_col_ind.flat[:]))
    fd_data = np.concatenate((fd_data, fd.flat[:]))

    current_row = row_ind[-1] + 1
    current_col = col_ind[-1] + 1

G = csr_matrix((fd_data, (row_ind, col_ind)), shape=(current_row, current_col))
save_sparse_npz('G.npz', G)
del G
G = load_sparse_npz('G.npz')


if False:
    G = G.toarray()
    import matplotlib.pyplot as plt
    plt.figure()
    plt.imshow(G) #.toarray())
    plt.show()
