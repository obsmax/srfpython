import numpy as np


# -------------------------------
def unpackmod96(string):
    """unpack 1D deptmodel depthdisp at mod96 format (see Herrmann's doc)
    """
    string = [line.strip() for line in string.split('\n')]
    if '' in string:
        string.remove('')
    # assert string[0].strip() == "MODEL.01"
    # title = string[1].strip()
    # isotropic = string[2].strip().upper() #== "ISOTROPIC"
    # kgs = string[3].strip().upper() #== "KGS"
    # flatearth = string[4].strip().upper() #== "FLAT EARTH"
    # oned = string[5].strip().upper() == "1-D"
    # cstvelo = string[6].strip().upper() == "CONSTANT VELOCITY"
    # assert string[7].strip() == "LINE08"
    # assert string[8].strip() == "LINE09"
    # assert string[9].strip() == "LINE10"
    # assert string[10].strip() == "LINE11"
    # header = string[11]

    assert string[11].startswith("H(KM)")
    nlayer = len(string) - 12
    #H, VP, VS, RHO, QP, QS, ETAP, ETAS, FREFP, FREFS = [np.empty(nlayer, float) for _ in xrange(10)]
    DAT = np.zeros((nlayer, 10), float)

    for n in xrange(nlayer):
        DAT[n, :] = np.asarray(string[12 + n].split(), float)

    H, VP, VS, RHO, QP, QS, ETAP, ETAS, FREFP, FREFS = [DAT[:, j] for j in xrange(10)]

    assert not H[-1]
    assert H[:-1].all()
    Z = np.concatenate(([0.], H[:-1].cumsum()))
    return nlayer, Z, H, VP, VS, RHO, QP, QS, ETAP, ETAS, FREFP, FREFS


# -------------------------------
def readmod96(filename):
    """read 1D deptmodel depthdisp at mod96 format (see Herrmann's doc)"""
    with open(filename, 'r') as fid:
        L = fid.readlines()
    return unpackmod96("".join(L))


# -------------------------------
def packmod96(Z, VP, VS, RHO, QP=None, QS=None, ETAP=None, ETAS=None, FREFP=None, FREFS=None):
    if QP is None: QP = np.zeros_like(VS)
    if QS is None: QS = np.zeros_like(VS)
    if ETAP is None: ETAP = np.zeros_like(VS)
    if ETAS is None: ETAS = np.zeros_like(VS)
    if FREFP is None: FREFP = np.ones_like(VS)
    if FREFS is None: FREFS = np.ones_like(VS)
    strout="""MODEL.01
MODELTITLE
ISOTROPIC
KGS
FLAT EARTH
1-D
CONSTANT VELOCITY
LINE08
LINE09
LINE10
LINE11
H(KM) VP(KM/S) VS(KM/S) RHO(GM/CC) QP QS ETAP ETAS FREFP FREFS
"""
    fmt = "%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n"

    H = np.zeros_like(Z)
    H[:-1] = Z[1:] - Z[:-1]
    H[-1] = 0.

    for tup in zip(H, VP, VS, RHO, QP, QS, ETAP, ETAS, FREFP, FREFS):
        strout += fmt % tup
    return strout
