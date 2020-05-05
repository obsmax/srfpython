from srfpython.depthdisp.dispcurves import surf96reader, surf96reader_from_surf96string, mklaws
import numpy as np
import os

"""
see theory.py
"""


def log_nofail(x):
    if np.isnan(x):
        return x
    elif x < 0.:
        return np.nan
    elif x == 0.:
        return -np.inf
    else:
        return np.log(x)


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
        # the default behavior is identity, see subclasses for advanced conversions
        d = values
        return d

    def inv(self, d):
        """converts a data array d into dispersion values"""
        # the default behavior is identity, see subclasses for advanced conversions
        values = d
        return values

    def inv_to_laws(self, d):
        values = self.inv(d)
        laws = mklaws(
            waves=self.waves, types=self.types,
            modes=self.modes, freqs=self.freqs,
            values=values, dvalues=self.dvalues)
        return laws


class Datacoder_log(Datacoder):
    def __init__(self, waves, types, modes, freqs, values, dvalues):
        Datacoder.__init__(self, waves, types, modes, freqs, values, dvalues)

    def target(self):
        dobs = np.log(self.values)
        CDinv = (self.dvalues / self.values) ** -2.
        return dobs, CDinv

    def __call__(self, values):
        assert len(values) == len(self.values)
        d = np.asarray(map(log_nofail, values), float)
        return d

    def inv(self, d):
        values = np.exp(d)
        return values


def makedatacoder(surf96filename, which=Datacoder_log):
    """
    :param surf96filename:
    :type surf96filename: str
    :param which: which datacoder to use
    :type which: type
    :return datacoder: a datacoder object initialized from the data found in surf96filename
    :rtype datacoder: DataCoder
    """
    if not isinstance(which, type):
        raise TypeError("which must be the class of the datacoder to use")

    if os.path.exists(surf96filename):
        s = surf96reader(surf96filename)
    else:
        s = surf96reader_from_surf96string(surf96filename)

    waves, types, modes, freqs, values, dvalues = s.wtmfvd()
    return which(waves, types, modes, freqs, values, dvalues)
