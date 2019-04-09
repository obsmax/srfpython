from srfpython.utils import discrete_time_primitive
from srfpython.depthdisp.dispcurves import surf96reader, groupbywtm
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
    def appendN(self, law, Ntimes=1):
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
    def appenddatN(self, f, v, Ntimes=1):
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
        coll, cax = ax.pcolormesh1(1. / self.f, self.v, H / I, **kwargs)
        ax.set_xscale('log')
        ax.set_yscale('log')
        return coll, cax

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
def dispstats(ds, percentiles=[0.16, 0.5, 0.84], Ndisp=100, weights=None, **mapkwargs):
    """
    compute statistics on a set of dispersion curves
    :param ds:
    :param percentiles:
    :param Ndisp:
    :param weights:
    :param mapkwargs:
    :return:
    """
    assert np.all([0 < p < 1 for p in percentiles])
    assert len(percentiles) == len(np.unique(percentiles))
    if weights is None:
        weights = np.ones(len(ds))
    else:
        assert len(weights) == len(ds)

    # initiate the depthpdfs
    waves, types, modes, freqs, values = ds[0]

    dpdfs = {}
    for w, t, m, f, _ in groupbywtm(waves, types, modes, freqs, values):
        dpdfs["%s%s%d" % (w, t, m)] = disppdf(f, np.logspace(np.log10(0.08), np.log10(4.0), Ndisp))

    for weight, (waves, types, modes, freqs, values) in zip(weights, ds):
        for w, t, m, f, v in groupbywtm(waves, types, modes, freqs, values):
            dpdfs["%s%s%d" % (w, t, m)].appenddatN(f, v, Ntimes=weight)

    for w, t, m, _, _ in groupbywtm(waves, types, modes, freqs, values):
        for p in percentiles:
            fpc, vpc = dpdfs["%s%s%d" % (w, t, m)].purcentile(p)
            wavespc = np.array([w for _ in xrange(len(fpc))], "|S1")
            typespc = np.array([t for _ in xrange(len(fpc))], "|S1")
            modespc = np.array([m for _ in xrange(len(fpc))], int)
            yield p, (wavespc, typespc, modespc, fpc, vpc)

