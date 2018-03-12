from srfpython.depthdisp.depthmodels import depthmodel1D, depthmodel
from srfpython.utils import discrete_time_primitive
import numpy as np
from srfpython.standalone.multipro8 import Job, StackAsync


# ----------------------------------------------------
class depthpdf(object):
    def __init__(self, z, v):
        assert np.all(z[1:] > z[:-1])
        assert np.all(v[1:] > v[:-1])
        self.z = z
        self.v = v
        self.H = np.zeros((len(z), len(v)), float)

    # ------------------------------------------------
    # def write(self, filename):
    #    assert filename.split('.')[-1] == "zpdf"
    #    pkl((self.z, self.v, self.H), filename)

    # ------------------------------------------------
    def append(self, m1D, **kwargs):
        # kwargs are passed to m1D.interp
        self.appendN(m1D, Ntimes=1, **kwargs)

    # ------------------------------------------------
    def appendN(self, m1D, Ntimes=1, **kwargs):
        """append the same model Ntimes times in the histogram"""

        def search_sorted_nearest(a, v, lena):
            n = np.searchsorted(a, v)
            if n == lena:
                return lena - 1
            elif n == 0:
                return 0
            else:
                vsup = a[n]
                vinf = a[n - 1]
                if abs(vsup - v) <= abs(v - vinf):
                    return n
                return n - 1

        zstairs, vstairs = m1D.stairs()
        zstairs[-1] = self.z.max() + 1.  # np.max([self.z.max(), 1.5 * zstairs[-2] + 1.0])
        lenz, lenv = len(self.z), len(self.v)
        # for i in xrange(len(zstairs)-1):
        # if i % 2: #horizontal step
        #     continue
        #     vmin, vmax = np.sort([vstairs[i], vstairs[i+1]])
        #     j0 = np.clip(np.searchsorted(self.v, vmin), 0, len(self.v) - 1)
        #     j1 = np.clip(np.searchsorted(self.v, vmax), 0, len(self.v))
        #     ii = np.clip(np.searchsorted(self.z, zstairs[i]), 0, len(self.z) - 1)
        #     self.H[ii,j0:j1] += float(Ntimes)
        # else: #vertical step
        #     zmin, zmax = zstairs[i], zstairs[i+1]
        #     i0 = search_sorted_nearest(self.z, zmin, lenz)#np.clip(np.searchsorted(self.z, zmin), 0, len(self.z) - 1)
        #     i1 = search_sorted_nearest(self.z, zmax, lenz) + 1#np.clip(np.searchsorted(self.z, zmax), 0, len(self.z))
        #     jj = search_sorted_nearest(self.v, vstairs[i], lenv) #np.clip(np.searchsorted(self.v, vstairs[i]), 0,len(self.v) - 1)
        #     self.H[i0:i1,jj] += float(Ntimes)

        for i in xrange(0, len(zstairs) - 1, 2):
            zmin, zmax = zstairs[i], zstairs[i + 1]
            i0 = search_sorted_nearest(self.z, zmin, lenz)  # np.clip(np.searchsorted(self.z, zmin), 0, len(self.z) - 1)
            i1 = search_sorted_nearest(self.z, zmax, lenz) + 1  # np.clip(np.searchsorted(self.z, zmax), 0, len(self.z))
            jj = search_sorted_nearest(self.v, vstairs[i],
                                       lenv)  # np.clip(np.searchsorted(self.v, vstairs[i]), 0,len(self.v) - 1)
            self.H[i0:i1, jj] += float(Ntimes)

    # ------------------------------------------------
    def show(self, ax, **kwargs):
        I = np.zeros(len(self.z))
        for i in xrange(len(self.z)):
            I[i] = discrete_time_primitive(self.v, self.H[i, :], area=True)
        I = discrete_time_primitive(self.z, I, area=True)

        H = np.ma.masked_where(self.H == 0., self.H)
        v_ = np.concatenate((self.v, [self.v[-1] + self.v[-1] - self.v[-2]]))
        z_ = np.concatenate((self.z, [self.z[-1] + self.z[-1] - self.z[-2]]))
        ax.pcolormesh1(v_, z_, H / I, **kwargs)
        if not ax.yaxis_inverted():
            ax.invert_yaxis()

    # ------------------------------------------------
    def purcentile(self, p):
        assert 0. <= p <= 1.
        P = np.zeros_like(self.z) * np.nan
        for i, z in enumerate(self.z):
            # x = self.v
            y = np.cumsum(self.H[i, :])
            ymin, ymax = y.min(), y.max()
            if ymax > ymin:
                y = (y - y.min()) / (y.max() - y.min())
                P[i] = np.interp(p, xp=y, fp=self.v)
            elif ymax == ymin:
                P[i] = self.v[0]
        assert not np.isnan(P).any()
        return self.z, P


# ------------------------------------------------
#class depthpdf_from_zpdffile(depthpdf):
#    def __init__(self, filename):
#        self.z, self.v, self.H = unpkl(filename)


# ------------------------------------------------
def dmstats(dms, percentiles=[0.16, 0.5, 0.84], Ndepth=100, Nvalue=100, weights=None):
    assert np.all([isinstance(dm, depthmodel) for dm in dms])
    assert np.all([0 < p < 1 for p in percentiles])
    assert len(percentiles) == len(np.unique(percentiles))
    if weights is None:
        weights = np.ones(len(dms))
    else:
        assert len(weights) == len(dms)

    zmax = -np.inf
    vsmin, vsmax = np.inf, -np.inf
    vpmin, vpmax = np.inf, -np.inf
    rhmin, rhmax = np.inf, -np.inf
    prmin, prmax = np.inf, -np.inf
    for dm in dms:
        zmax = np.max([zmax, 1.1 * dm.vs.z[-1]])
        vsmin = np.min([vsmin, dm.vs.values.min()])
        vsmax = np.max([vsmax, dm.vs.values.max()])
        vpmin = np.min([vpmin, dm.vp.values.min()])
        vpmax = np.max([vpmax, dm.vp.values.max()])
        rhmin = np.min([rhmin, dm.rh.values.min()])
        rhmax = np.max([rhmax, dm.rh.values.max()])
        prmin = np.min([prmin, dm.pr().values.min()])
        prmax = np.max([prmax, dm.pr().values.max()])

    zbins = np.linspace(0., zmax, Ndepth)
    vspdf = depthpdf(z=zbins, v=np.linspace(vsmin, vsmax, Nvalue))
    vppdf = depthpdf(z=zbins, v=np.linspace(vpmin, vpmax, Nvalue))
    rhpdf = depthpdf(z=zbins, v=np.linspace(rhmin, rhmax, Nvalue))
    prpdf = depthpdf(z=zbins, v=np.linspace(prmin, prmax, Nvalue))

    for dm, weight in zip(dms, weights):
        vspdf.appendN(dm.vs, Ntimes=weight)
        vppdf.appendN(dm.vp, Ntimes=weight)
        rhpdf.appendN(dm.rh, Ntimes=weight)
        prpdf.appendN(dm.pr(), Ntimes=weight)

    for p in percentiles:
        zpc, vspc = vspdf.purcentile(p)
        _, vppc = vppdf.purcentile(p)
        _, rhpc = rhpdf.purcentile(p)
        _, prpc = prpdf.purcentile(p)
        vppc = depthmodel1D(zpc, vppc)  # .simplify()
        vspc = depthmodel1D(zpc, vspc)  # .simplify()
        rhpc = depthmodel1D(zpc, rhpc)  # .simplify()
        prpc = depthmodel1D(zpc, prpc)  # .simplify()

        # dmpc = depthmodel(\).simplify()
        # yield p, dmpc #, (zmed, vpmed, vsmed, rhmed, prmed)
        yield p, (vppc, vspc, rhpc, prpc)


# ------------------------------------------------
class UserStacker(object):
    def __init__(self, zmax, Ndepth, vpmin, vpmax, vsmin, vsmax, rhmin, rhmax, prmin, prmax, Nvalue):
        zbins = np.linspace(0., zmax, Ndepth)
        self.vspdf = depthpdf(z=zbins, v=np.unique(np.linspace(vsmin, vsmax, Nvalue)))
        self.vppdf = depthpdf(z=zbins, v=np.unique(np.linspace(vpmin, vpmax, Nvalue)))
        self.rhpdf = depthpdf(z=zbins, v=np.unique(np.linspace(rhmin, rhmax, Nvalue)))
        self.prpdf = depthpdf(z=zbins, v=np.unique(np.linspace(prmin, prmax, Nvalue)))

    def __call__(self, worker, weight, dm):
        self.vspdf.appendN(dm.vs, Ntimes=weight)
        self.vppdf.appendN(dm.vp, Ntimes=weight)
        self.prpdf.appendN(dm.pr(), Ntimes=weight)
        self.rhpdf.appendN(dm.rh, Ntimes=weight)
        return self, worker.name

    def __iadd__(self, other):
        "a method to merge stackers once back to the serial section"
        assert isinstance(other, UserStacker)
        self.vppdf.H += other.vppdf.H
        self.vspdf.H += other.vspdf.H
        self.prpdf.H += other.prpdf.H
        self.rhpdf.H += other.rhpdf.H
        return self


def dmstats1(dms, percentiles=[0.16, 0.5, 0.84], Ndepth=100, Nvalue=100, weights=None, **mapkwargs):
    assert np.all([isinstance(dm, depthmodel) for dm in dms])
    assert np.all([0 < p < 1 for p in percentiles])
    assert len(percentiles) == len(np.unique(percentiles))
    if weights is None:
        weights = np.ones(len(dms))
    else:
        assert len(weights) == len(dms)

    zmax = -np.inf
    vsmin, vsmax = np.inf, -np.inf
    vpmin, vpmax = np.inf, -np.inf
    rhmin, rhmax = np.inf, -np.inf
    prmin, prmax = np.inf, -np.inf
    for dm in dms[:10]:
        zmax = np.max([zmax, 1.1 * dm.vs.z[-1]])
        vsmin = np.min([vsmin, dm.vs.values.min()])
        vsmax = np.max([vsmax, dm.vs.values.max()])
        vpmin = np.min([vpmin, dm.vp.values.min()])
        vpmax = np.max([vpmax, dm.vp.values.max()])
        rhmin = np.min([rhmin, dm.rh.values.min()])
        rhmax = np.max([rhmax, dm.rh.values.max()])
        prmin = np.min([prmin, dm.pr().values.min()])
        prmax = np.max([prmax, dm.pr().values.max()])
    if vsmin == vsmax:
        vsmin *= 0.99
        vsmax *= 1.01
    if vpmin == vpmax:
        vpmin *= 0.99
        vpmax *= 1.01
    if rhmin == rhmax:
        rhmin *= 0.99
        rhmax *= 1.01
    if prmin == prmax:
        prmin *= 0.99
        prmax *= 1.01

    # ----------------------
    def JobGen():
        for weight, dm in zip(weights, dms):
            yield Job(weight, dm)

    # ----------------------
    s0 = UserStacker(zmax, Ndepth, vpmin, vpmax, vsmin, vsmax, rhmin, rhmax, prmin, prmax,
                     Nvalue)  # the stacker to be reproduced in each independent workspace (deep copy)
    with StackAsync(s0, JobGen(), Verbose=False, **mapkwargs) as sa:
        S = UserStacker(zmax, Ndepth, vpmin, vpmax, vsmin, vsmax, rhmin, rhmax, prmin, prmax,
                        Nvalue)  # initiate the final stack
        for jobids, (s, wname), Tgen, Tpro in sa:  # receive partial stacks
            sa.communicate("Stacker %s stacked %6d jobs in %.6fs" % (wname, len(jobids), Tgen + Tpro))
            S += s  # merge all partial stacks using method __iadd__

    for p in percentiles:
        zpc, vspc = S.vspdf.purcentile(p)
        _, vppc = S.vppdf.purcentile(p)
        _, rhpc = S.rhpdf.purcentile(p)
        _, prpc = S.prpdf.purcentile(p)

        vppc = depthmodel1D(zpc, vppc)  # .simplify()
        vspc = depthmodel1D(zpc, vspc)  # .simplify()
        rhpc = depthmodel1D(zpc, rhpc)  # .simplify()
        prpc = depthmodel1D(zpc, prpc)  # .simplify()

        # dmpc = depthmodel(\).simplify()
        # yield p, dmpc #, (zmed, vpmed, vsmed, rhmed, prmed)
        yield p, (vppc, vspc, rhpc, prpc)
