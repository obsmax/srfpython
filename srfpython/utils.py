from __future__ import print_function
from numpy import log10
import time
import signal
import numpy as np


class TimeOutError(Exception):
    pass


class Timeout():

    def __init__(self, sec):
        self.sec = sec

    def __enter__(self):
        signal.signal(signal.SIGALRM, self.raise_timeout)
        signal.alarm(self.sec)

    def __exit__(self, *args):
        signal.alarm(0)

    def raise_timeout(self, *args):
        raise TimeOutError()


class Timer(object):

    def __init__(self, title):
        self.title = title

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args, **kwargs):
        print ("elapsed time %s : %fs" % (self.title, time.time() - self.start))


def firstfalse(I):
    if not I[0] : return 0
    II = I[1:] != I[:-1]
    r = np.arange(len(I) - 1)
    return r[II][0] + 1


def munique(*Xs):
    assert np.all(len(X) == len(Xs[0]) for X in Xs)
    L = []
    for tup in zip(*Xs):
        if tup not in L:
            L.append(tup)
    return tuple([np.array(w) for w in zip(*L)])


def freqspace(freqmin, freqmax, nfreq, scale="flin"):
    if "lin" in scale.lower():
        return np.linspace(freqmin, freqmax, nfreq)
    elif "log" in scale.lower():
        return np.logspace(log10(freqmin), log10(freqmax), nfreq)
    else:
        raise ValueError('%s not understood' % scale)


def minmax(X):
    if hasattr(X, "min"): #numpy arrays
        return X.min(), X.max()
    else:
        return min(X), max(X)


def discrete_time_primitive(x, u, area=False):
    if len(x) - len(u):
        raise Exception('shape mismatch')
    if (x[1:] <= x[:-1]).any():
        raise Exception('x must be sorted asc')

    v = np.zeros(len(u), u.dtype)
    v[1:] = 0.5 * (u[1:] + u[:-1]) * (x[1:] - x[:-1])
    if area:
        return np.sum(v)  # area below the curve
    else:
        return v.cumsum()  # primitive function to be associated to x


def cosTaperwidth(data, sampling_rate, width):
    tap = np.ones_like(data)
    Nwidth = int(np.round(width * sampling_rate))
    if not Nwidth : return tap

    t = np.arange(len(data)) / sampling_rate
    ttap = 0.5 * (np.sin(np.pi * t[:Nwidth] / float(width) / 1.0 + np.pi / 2.0) + 1.)
    tap[:Nwidth] *= ttap[::-1]
    tap[-Nwidth:] *= ttap
    return tap
