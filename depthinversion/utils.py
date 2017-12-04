from pickle import load, dump
import numpy as np
import time



class Timer(object):
    def __init__(self, title):
        self.title = title
    def __enter__(self):
        self.start = time.time()
        return self
    def __exit__(self, *args, **kwargs):
        print ("elapsed time %s : %fs" % (self.title, time.time() - self.start))

# -------------------------
def pkl(obj, filename, protocol = None):
    with open(filename, 'wb') as fid:
        dump(obj, fid, protocol = protocol)

# -------------------------
def unpkl(filename):
    with open(filename, 'rb') as fid:
        obj = load(fid)
    return obj

# -------------------------
def cosTaperwidth(data, sampling_rate, width):
    tap = np.ones_like(data)
    Nwidth = int(np.round(width * sampling_rate))
    if not Nwidth : return tap

    t = np.arange(len(data)) / sampling_rate
    ttap = 0.5 * (np.sin(np.pi * t[:Nwidth] / float(width) / 1.0 + np.pi / 2.0) + 1.)
    tap[:Nwidth]  *= ttap[::-1]
    tap[-Nwidth:] *= ttap
    return tap

# -------------------------
def discrete_time_primitive(x, u, area=False):
    if len(x) - len(u):
        raise Exception('shape mismatch')
    if not (x[1:] > x[:-1]).all():
        raise Exception('x must be sorted asc')

    v = np.zeros(len(u), u.dtype)
    v[1:] = 0.5 * (u[1:] + u[:-1]) * (x[1:] - x[:-1])
    if area:
        return np.sum(v) #area below the curve
    else:
        return v.cumsum() #primitive function to be associated to x
# -------------------------
def munique(*Xs):
    """ multiple unique
    same as np.unique but for tuples
    e.g.
        x = [1, 2, 3, 1, 1]
        y = [2, 3, 2, 2, 3]
    munique(x, y)
    returns np.array([1, 2, 3, 1])
            np.array([2, 3, 2, 3])

    the pair (1, 2) has been removed since it appeared twice
    the pairs (3, 2) and (2, 3) are kept
    """
    assert np.all(len(X) == len(Xs[0]) for X in Xs[1:])
    L = []
    for tup in zip(*Xs):
        if tup not in L:
            L.append(tup)
    return tuple([np.array(w) for w in zip(*L)])
# -------------------------
def freqspace(freqmin, freqmax, nfreq, scale="flin"):
    if "lin" in scale.lower():
        return np.linspace(freqmin, freqmax, nfreq)
    elif "log" in scale.lower():
        return np.logspace(log10(freqmin), log10(freqmax), nfreq)
    else: raise ValueError('%s not understood' % scale)

# -------------------------
def histogram2d(xflat, yflat, xbins, ybins):
    H, X, Y = np.histogram2d(x = xflat, y = yflat, bins=(xbins, ybins), normed=True)
    #H is in bin^-1 * xunit ^ -1 * yunit ^ -1
    #H is a true pdf, indeed (H * dX * dY).sum() = 1.0
    #however H can exceed 1. locally if the bin area is lower than 1.0 !!!!!
    H[np.isnan(H)] = 0.
    return X, Y, H.T