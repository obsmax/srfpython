import numpy as np


def unpacksurf96(string):
    """unpack dispersion curves at surf96 format (see Herrmann's doc)"""
    string = [line.strip() for line in string.split('\n')]
    if "" in string:
        string.remove('')
    npoints = len(string)

    datatypes = ['|S1', '|S1', '|S1', int, float, float, float]
    WAVE, TYPE, FLAG, MODE, PERIOD, VALUE, DVALUE = [np.empty(npoints, dtype=d) for d in datatypes]
    NLC, NLU, NRC, NRU = 0, 0, 0, 0
    for n in range(npoints):
        l = string[n].split()
        WAVE[n], TYPE[n], FLAG[n] = np.asarray(l[1:4], "|S1")
        MODE[n] = int(l[4])
        PERIOD[n], VALUE[n], DVALUE[n] = np.asarray(l[5:], float)
        if WAVE[n] == "L":
            if TYPE[n] == "C":
                NLC += 1
            elif TYPE[n] == "U":
                NLU += 1
        elif WAVE[n] == "R":
            if TYPE[n] == "C":
                NRC += 1
            elif TYPE[n] == "U":
                NRU += 1
        else:
            raise Exception('')
    return WAVE, TYPE, FLAG, MODE, PERIOD, VALUE, DVALUE, NLC, NLU, NRC, NRU


def readsurf96(filename):
    """read dispersion depthdisp at surf96 format"""
    with open(filename, 'r') as fid:
        L = fid.readlines()
    return unpacksurf96("".join(L))


def packsurf96(wave, type, flag, mode, period, value, dvalue):
    fmt = "SURF96 {wave} {type} {flag} {mode} {period} {value} {dvalue}"
    assert np.all([len(_) == len(wave) for _ in [type, flag, mode, period, value, dvalue]])
    mode = np.asarray(mode, int)

    return "\n".join([fmt.format(
        wave=w, type=t, mode=m, period=p,
        value=v, dvalue=d, flag=f)
        for w, t, f, m, p, v, d in
            zip(wave, type, flag, mode, period, value, dvalue)])
