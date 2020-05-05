import numpy as np
from scipy.sparse import diags, csr_matrix, csc_matrix, \
    save_npz as save_sparse_npz, load_npz as load_sparse_npz
from scipy.sparse.linalg import inv as sparseinv
from srfpython.HerrMet.paramfile import load_paramfile
from srfpython.depthdisp.depthmodels import depthmodel_from_arrays
from srfpython.HerrMet.datacoders import makedatacoder, Datacoder_log
from srfpython.HerrMet.theory import Theory
from srfpython.standalone.display import *


CM = load_sparse_npz('CM.npz')
G = load_sparse_npz('G.npz')
CD = load_sparse_npz('CD.npz')
Dobs = np.load('Dobs.npy')
Mprior = np.load('Mprior.npy')
paramterizer_strings = np.load('paramterizer_strings.npy')
datacoder_strings = np.load('datacoder_strings.npy')
parameterizers = [load_paramfile(paramterizer_string, verbose=False)[0]
                  for paramterizer_string in paramterizer_strings]
datacoders = [makedatacoder(datacoder_string, Datacoder_log)
                  for datacoder_string in datacoder_strings]
theorys = [Theory(parameterizer=parameterizer, datacoder=datacoder)
           for parameterizer, datacoder in zip(parameterizers, datacoders)]

def g(m):
    ms = m.reshape((len(m) // parameterizers[0].NLAYER, parameterizers[0].NLAYER))
    ds = []
    for n in range(ms.shape[0]):
        ds.append(theorys[n](ms[n, ...]))
    return np.concatenate(ds)


CMGT = CM * G.T
print(CMGT.shape)
# print(np.sum(CMGT == 0) / float(CMGT.size))
S = CD + G * CMGT
Sinv = sparseinv(S.tocsc())

M0 = Mprior

MS = [M0]
MI = M0
for _ in range(2):
    XI = Dobs - g(MI) + G * (MI - Mprior)
    MI = M0 + CMGT * Sinv * XI
    MS.append(MI.copy())


def MtoDM(M3):
    print(type(M3))
    print(M3.shape)
    print(parameterizers[0].NLAYER)
    ms = M3.reshape((len(M3) // parameterizers[0].NLAYER, parameterizers[0].NLAYER))
    dms = []
    for n in range(ms.shape[0]):
        ztop, vp, vs, rh = parameterizers[n].inv(ms[n, ...])
        print (n, vs)
        dms.append(depthmodel_from_arrays(ztop, vp, vs, rh))
    return dms


dms0 = MtoDM(MS[0])
dmslast = MtoDM(MS[-1])

def showdms(ax, dms, **kwargs):
    vs = []
    for n, dm in enumerate(dms):
        vs.append(dm.vs.values)
    x = np.arange(len(vs)+1) - 0.5
    y = np.hstack((dm.vs.z, dm.vs.z[-1] + 0.1))
    plt.colorbar(ax.pcolormesh(
        x,
        y,
        np.array(vs).T,
        **kwargs), ax=ax)

def showdiff(ax, dms1, dms2, **kwargs):
    vsdiff = []
    for n, (dm1, dm2) in enumerate(zip(dms1, dms2)):
        vsdiff.append(dm1.vs.values - dm2.vs.values)
    x = np.arange(len(vsdiff)+1) - 0.5
    y = np.hstack((dm1.vs.z, dm1.vs.z[-1] + 0.1))
    plt.colorbar(ax.pcolormesh(
        x,
        y,
        np.array(vsdiff).T,
        **kwargs), ax=ax)


ax0 = plt.subplot(311)
ax3 = plt.subplot(312, sharex=gca(), sharey=gca())
axdiff = plt.subplot(313, sharex=gca(), sharey=gca())
showdms(ax0,
        dms0,
        vmin=0.5, vmax=3.5, cmap=plt.get_cmap('nipy_spectral_r'))
showdms(ax3,
        dmslast,
        vmin=0.5, vmax=3.5, cmap=plt.get_cmap('nipy_spectral_r'))

showdiff(axdiff,
         dms0, dmslast,
         cmap=plt.get_cmap('RdBu'))

ax0.invert_yaxis()
showme()





