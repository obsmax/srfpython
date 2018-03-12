import numpy as np


# -------------------------------
def unpackmod96(string):
    """unpack 1D deptmodel depthdisp at mod96 format (see Herrmann's doc)
    """
    string = [line.strip() for line in string.split('\n')]

    # remove blank lines
    # string.remove('') #fucking not working all the time
    I = [_ != '' for _ in string]
    string = [string[i] for i in np.arange(len(string))[I]]

    # make sure the header is correctly formated
    assert string[0] == "MODEL.01"
    # title = string[1].strip()
    assert string[2].upper() == "ISOTROPIC"
    assert string[3].upper() == "KGS"
    assert string[4].upper() == "FLAT EARTH"
    assert string[5].upper() == "1-D"
    assert string[6].upper() == "CONSTANT VELOCITY"
    assert string[7] == "LINE08"
    assert string[8] == "LINE09"
    assert string[9] == "LINE10"
    assert string[10] == "LINE11"
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
