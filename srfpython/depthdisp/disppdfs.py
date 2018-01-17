from srfpython.utils import discrete_time_primitive
import numpy as np


# ----------------------------------------------------
class disppdf(object):
    #2D histogram for dipsersion curves, see also depthmodels.depthpdf, depthmodels.depthpdf2disppdf
    def __init__(self, f, v):
        assert np.all(f[1:] > f[:-1])
        assert np.all(v[1:] > v[:-1])
        self.f = f
        self.v = v
        self.H = np.zeros((len(v), len(f)), float)

    # ------------------------------------------------
    # def write(self, filename):
    #    assert filename.split('.')[-1] == "dpdf"
    #    pkl((self.f, self.v, self.H), filename)

    # ------------------------------------------------
    def append(self, law):
        self.appendN(law, Ntimes = 1)

    # ------------------------------------------------
    def appendN(self, law, Ntimes=1, **kwargs):
        """append the same model Ntimes times in the histogram"""
        v = law(self.f)
        for i in xrange(len(self.f)):
            if np.isnan(v[i]):continue
            j = np.clip(np.searchsorted(self.v, v[i]), 0, len(self.v) - 1)
            self.H[j, i] += float(Ntimes)
#            if 0 <= j < len(self.v):
#                self.H[j, i] += float(Ntimes)

    # ------------------------------------------------
    def appenddat(self, f, v):
        self.appenddatN(f, v, Ntimes = 1)

    # ------------------------------------------------
    def appenddatN(self, f, v, Ntimes=1, **kwargs):
        if (f != self.f).any():
            v = np.interp(self.f, xp = f, fp = v, left = np.nan, right = np.nan)
        for i in xrange(len(self.f)):
            if np.isnan(v[i]):continue
            j = np.clip(np.searchsorted(self.v, v[i]), 0, len(self.v) - 1)
            self.H[j, i] += float(Ntimes)

    # ------------------------------------------------
    def show(self, ax, **kwargs):
        I = np.zeros(len(self.f))
        for i in xrange(len(self.f)):
            I[i] = discrete_time_primitive(self.v, self.H[:, i], area = True)
        I = discrete_time_primitive(self.f, I, area = True)
        H = np.ma.masked_where(self.H == 0., self.H)
        ax.pcolormesh1(1. / self.f, self.v, H / I, **kwargs)
        ax.set_xscale('log')
        ax.set_yscale('log')

    # ------------------------------------------------
    def purcentile(self, p):
        assert 0. <= p <= 1.
        P = np.zeros_like(self.f) * np.nan
        for j, f in enumerate(self.f):
            x = self.v
            y = np.cumsum(self.H[:, j])
            ymin, ymax = y.min(), y.max()
            if ymax > ymin:
                y = (y - y.min()) / (y.max() - y.min())
                P[j] = np.interp(p, xp = y, fp = x)
        I = ~np.isnan(P)
        return self.f[I], P[I]


# ------------------------------------------------
# class disppdf_from_zpdffile(disppdf):
#    def __init__(self, filename):
#        self.f, self.v, self.H = unpkl(filename)
