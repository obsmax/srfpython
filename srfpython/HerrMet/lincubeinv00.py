from __future__ import print_function
from builtins import input
import sys, glob, os
import numpy as np
import matplotlib.pyplot as plt
from srfpython.HerrMet.nodefile import NodeFileString, NodeFile
from srfpython.HerrMet.theory import Theory
from srfpython.HerrMet.paramfile import load_paramfile
from srfpython.HerrMet.datacoders import makedatacoder, Datacoder_log
from srfpython.standalone.asciifile import AsciiFile
from srfpython.HerrMet.files import HERRMETEXTRACTPDFMODELFILE, HERRMETTARGETFILE, ROOTNAME
from srfpython.depthdisp.depthmodels import depthmodel_from_mod96, brocher2005, depthmodel_from_arrays, depthmodel1D, depthspace
from scipy.sparse import diags, csr_matrix, csc_matrix, save_npz as save_sparse_npz, load_npz as load_sparse_npz

"""
set the parameterization, 
load the prior model and uncertainty
load the data to fit and uncertainties
"""

nodefile = "../../tutorials/02_cube_inversion_example/nodes.txt"
# resample the aposteri median at ztop defined below
ztop = np.linspace(0., 3.0, 50)
horizontal_smoohting_distance = 50.0 * np.ones_like(ztop)  # one smoothing distance per layer
vertical_smoothing_distance = 2.0
kkk = 50.

# node file
nf = NodeFile(nodefile)
nf.fill_extraction_files()

zmid = ztop + np.hstack(
    (0.5 * (ztop[1:] - ztop[:-1]),
     0.5 * np.median(ztop[1:] - ztop[:-1])))

def haversine(loni, lati, lonj, latj, R=6371.):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    source https://stackoverflow.com/questions/15736995/how-can-i-quickly-estimate-the-distance-between-two-latitude-longitude-points
    lone, late : coordinates of the first point, in degrees
    lons, lats : coordinates of the second point, in degrees
    :return distance in km

    consistent with Ll2DGAB
    """
    # convert decimal degrees to radians
    q = np.pi / 180.

    # haversine formula
    dlon = (loni - lonj)
    dlat = (lati - latj)
    a = np.sin(q * dlat / 2.) ** 2. + np.cos(q * latj) * np.cos(q * lati) * np.sin(q * dlon / 2.) ** 2.
    c = 2. * np.arcsin(np.sqrt(a))
    km = R * c
    return km

# ============================= G, CDinv, Dobs, Dcalc
CDinv_diag_data = np.array([], float)
Dobs = np.array([], float)
Mprior = np.array([], float)
parameterizer_strings = []
datacoder_strings = []
for nnode, (node, lon, lat, targetfile, paramfile, medianfile, p16file, p84file) in enumerate(nf):

    print(node, lon, lat)
    print("    ", targetfile)
    print("    ", paramfile)
    print("    ", medianfile)
    print("    ", p16file)
    print("    ", p84file)

    # load the extracted file
    dm = depthmodel_from_mod96(medianfile)

    # resample to ztop
    vs = dm.vs.interp(zmid)
    # vs = np.linspace(0.8, 3.5, len(vs))
    # plt.figure()
    # plt.plot(vs, zmid)
    # plt.gca().invert_yaxis()
    # plt.show()
    # exit()

    # ============================= G
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
        # force VINF=VSUP => means lock the depth of the interfaces in the theory operator
        parameter_string += "-Z{} {} {}\n".format(i, -ztop[i], -ztop[i])  # add locked depth interfaces

    for i in range(len(ztop)):
        # SET VINF < VS extracted from pointwise inv < VSUP
        # such as parameterizer.MMEAN corresponds to the extracted vs
        parameter_string += "VS{} {} {}\n".format(i, vs[i]-0.01, vs[i]+0.01)

    # initiate the parameterizer and datacoder
    parameterizer_strings.append(parameter_string)
    with open(targetfile, 'r') as fid:
        datacoder_string = "".join(fid.readlines())
        datacoder_strings.append(datacoder_string)
    parameterizer, _ = load_paramfile(parameter_string, verbose=False)
    datacoder = makedatacoder(targetfile, which=Datacoder_log)

    # ============================= Mprior, CM
    mnode = parameterizer.MMEAN
    Mprior = np.concatenate((Mprior, mnode))

    # ============================= CDinv and Dobs
    dobs_current, CDinv_diag_current = datacoder.target()
    CDinv_diag_data = np.concatenate((CDinv_diag_data, CDinv_diag_current))
    Dobs = np.concatenate((Dobs, dobs_current))

# ============================= CM
CM_row_ind = np.array([], int)
CM_col_ind = np.array([], int)
CM_data = np.array([], float)
for nnode in range(len(nf)):
    # find the posterior unc at each depth on vs from the pointwise inversion

    vs84_n = depthmodel_from_mod96(nf.p84files[nnode]).vs
    vs16_n = depthmodel_from_mod96(nf.p16files[nnode]).vs
    vs84_n = depthmodel1D(ztop, vs84_n.interp(zmid))
    vs16_n = depthmodel1D(ztop, vs16_n.interp(zmid))
    vs_unc_n = 0.5 * (vs84_n.values - vs16_n.values) * kkk

    # determine the vertical correlation coeff in cell n
    if vertical_smoothing_distance > 0:
        rhonn = np.exp(-0.5 * ((ztop - ztop[:, np.newaxis]) / vertical_smoothing_distance) ** 2.)
        covnn = vs_unc_n[:, np.newaxis] * vs_unc_n * rhonn

        row_ind, col_ind = np.meshgrid(range(len(ztop)), range(len(ztop)))
        CM_row_ind = np.hstack((CM_row_ind, nnode * len(ztop) + row_ind.flat[:]))
        CM_col_ind = np.hstack((CM_col_ind, nnode * len(ztop) + col_ind.flat[:]))
        CM_data = np.hstack((CM_data, covnn.flat[:]))
    else:
        covnn_diags = vs_unc_n ** 2.0
        CM_row_ind = np.hstack((CM_row_ind, nnode * len(ztop) + np.arange(len(ztop))))
        CM_col_ind = np.hstack((CM_col_ind, nnode * len(ztop) + np.arange(len(ztop))))
        CM_data = np.hstack((CM_data, covnn_diags.flat[:]))

    for mnode in range(nnode + 1, len(nf)):
        lonn = nf.lons[nnode]
        latn = nf.lats[nnode]
        lonm = nf.lons[mnode]
        latm = nf.lats[mnode]

        vs84_m = depthmodel_from_mod96(nf.p84files[mnode]).vs
        vs16_m = depthmodel_from_mod96(nf.p16files[mnode]).vs
        vs84_m = depthmodel1D(ztop, vs84_m.interp(zmid))
        vs16_m = depthmodel1D(ztop, vs16_m.interp(zmid))
        vs_unc_m = 0.5 * (vs84_m.values - vs16_m.values) * kkk

        dnm = haversine(lonn, latn, lonm, latm)
        rhonm_diags = np.exp(-0.5 * (dnm / horizontal_smoohting_distance) ** 2.)
        covnm_diags = vs_unc_n * vs_unc_m * rhonm_diags

        CM_row_ind = np.hstack((CM_row_ind, nnode * len(ztop) + np.arange(len(ztop))))
        CM_col_ind = np.hstack((CM_col_ind, mnode * len(ztop) + np.arange(len(ztop))))
        CM_data = np.hstack((CM_data, covnm_diags.flat[:]))

        CM_row_ind = np.hstack((CM_row_ind, mnode * len(ztop) + np.arange(len(ztop))))
        CM_col_ind = np.hstack((CM_col_ind, nnode * len(ztop) + np.arange(len(ztop))))
        CM_data = np.hstack((CM_data, covnm_diags.flat[:]))

CM = csc_matrix((CM_data, (CM_row_ind, CM_col_ind)), shape=(len(nf) * len(ztop), len(nf) * len(ztop)))
if 0:
    assert input('sure ({})?)'.format(CM.shape)) == "y"
    CM = CM.toarray()
    assert (CM.T == CM).all()
    plt.figure()
    plt.imshow(CM)
    plt.show()
    exit()


# ============================= save matrixes to disk
CD = diags(CDinv_diag_data ** -1.0, offsets=0, format="csr", dtype=float)
CDinv = diags(CDinv_diag_data, offsets=0, format="csr", dtype=float)
parameterizer_strings = np.asarray(parameterizer_strings, str)
datacoder_strings = np.asarray(datacoder_strings, str)
# np.save('CDinv_diag.npy', CDinv_diag, allow_pickle=False)
save_sparse_npz('CD.npz', CDinv)
save_sparse_npz('CDinv.npz', CDinv)
save_sparse_npz('CM.npz', CM)
np.save('Dobs.npy', Dobs, allow_pickle=False)
np.save('Mprior.npy', Mprior, allow_pickle=False)
np.save('paramterizer_strings.npy', parameterizer_strings, allow_pickle=False)
np.save('datacoder_strings.npy', datacoder_strings, allow_pickle=False)
