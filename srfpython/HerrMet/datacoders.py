from srfpython.depthdisp.dispcurves import surf96reader, surf96reader_from_surf96string
import numpy as np
import os

"""
see theory.py
"""


class Datacoder(object):
    def __init__(self, waves, types, modes, freqs, values, dvalues):
        """init with the target data and uncertainty"""

        assert np.all(~np.isnan(values))

        self.npoints = len(waves)
        assert len(types) == len(modes) == \
               len(freqs) == len(values) == \
               len(dvalues) == self.npoints

        self.waves = waves      # parameters for forward problem
        self.types = types      # parameters for forward problem
        self.modes = modes      # parameters for forward problem
        self.freqs = freqs      # parameters for forward problem
        self.values = values    # target dispersion
        self.dvalues = dvalues  # target dispersion

    def target(self):
        dobs = self.values  # target data array
        CDinv = self.dvalues ** -2.  # target data covariance (inverted, diagonal terms)
        return dobs, CDinv

    def __call__(self, values):
        """converts dispersion values into a data array d"""
        # assert len(values) == self.npoints
        return values  # the default behavior is identity, see subclasses for advanced conversions

    def inv(self, d):
        """converts a data array d into dispersion values"""
        return d  # the default behavior is identity, see subclasses for advanced conversions


def log_nofail(x):
    if np.isnan(x):
        return x
    elif x < 0.:
        return np.nan
    elif x == 0.:
        return -np.inf
    else:
        return np.log(x)


class Datacoder_log(Datacoder):
    def __init__(self, waves, types, modes, freqs, values, dvalues):
        Datacoder.__init__(self, waves, types, modes, freqs, values, dvalues)

    def target(self):
        dobs   = np.log(self.values)
        CDinv  = (self.dvalues / self.values) ** -2.
        return dobs, CDinv

    def __call__(self, values):
        assert len(values) == len(self.values)
        return np.asarray(map(log_nofail, values), float)
        # d = np.log(values)
        # return d

    def inv(self, d):
        values = np.exp(d)
        return values


def makedatacoder(s96, which=Datacoder_log):
    if os.path.exists(s96):
        s = surf96reader(s96)
    else:
        s = surf96reader_from_surf96string(s96)

    waves, types, modes, freqs, values, dvalues = s.wtmfvd()
    return which(waves, types, modes, freqs, values, dvalues)
