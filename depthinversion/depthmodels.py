from ..processingtools.multipro7 import Job, StackAsync
from .utils import cosTaperwidth, pkl, unpkl, discrete_time_primitive
from scipy.fftpack import fft, ifft, fftfreq
import numpy as np
import copy, os
# from scipy.misc import comb


#----------------------------------------------------
def gardner74(vp):
    """Gardner's density law for sedimentary rocks, 1974
       vp between 1.5 and 6.1 km/s
    """
    return 1.74 * vp ** 0.25 #density in g/cm3
#----------------------------------------------------
def depthspace(zbot, n):
    """
    zbot  = top of the half space
    n     = desired number of layers
    """
    z = np.logspace(np.log10(0.1), np.log10(3.), n)
    z = zbot * (z - z.min()) / (z.max() - z.min()) 
    return z
#----------------------------------------------------
def randab(r, a, b):
    #assert 0. <= r <= 1.
    return r * (b - a) + a
#----------------------------------------------------
class depthmodel1D(object):
    """self.z = ztop !!!
       can handle irregular layers
    """
    #----------------------------------------------------
    def __init__(self, ztop, values, interpmethod = "stairs"):
        #last value of v stand for half space
        assert len(ztop) == len(values)
        assert np.all(ztop[1:] > ztop[:-1])
        assert ztop[0] == 0.
        self.values, self.z = np.array(values), np.array(ztop)
        self.interpmethod = interpmethod #to be used by default

    #----------------------------------------------------
    def __len__(self):
        return len(self.z)#number of layers including the lower half space

    #----------------------------------------------------
    def ztop(self): return self.z

    #----------------------------------------------------
    def thickness(self): return np.concatenate((self.z[1:] - self.z[:-1], [np.inf]))

    #----------------------------------------------------
    def zbot(self):  return self.z + self.thickness()

    #----------------------------------------------------
    def zmid(self):  return self.z + 0.5 * self.thickness()

    #----------------------------------------------------
    def stairs(self):
        """
        :return:
            a stairs like version of the model for display or interpolation
        """
        N = len(self)
        Z = self.zbot()
        I = np.sort( np.concatenate( (range(N-1), range(N)) ) )
        Z = np.concatenate( ([0.], Z[I]) )
        I = np.sort( np.concatenate( (range(N), range(N)) ) )
        return Z, self.values[I]

    #----------------------------------------------------
    def show(self, ax, marker = "-", **kwargs): #halfspacemarker="-", 

        zz, vv = self.stairs()
        zz[-1] = np.max([1.5 * self.z[-1], 3.0]) #cannot plot np.inf...

        hdl = ax.plot(vv, zz, marker, **kwargs)[0]
        #hdl = ax.plot(vv[:-1], zz[:-1], marker, **kwargs)[0]
        #kwargs["linestyle"] = halfspacemarker
        #ax.plot(vv[-2:], zz[-2:], marker, **kwargs)

        if not ax.yaxis_inverted():
            ax.invert_yaxis()
        return hdl

    #----------------------------------------------------
    def interp(self, z, interpmethod = None):
        """
        :param z: depth at which to interpolate the model
        :param interpmethod:
            method "stairs" means simply reading the "staired" model at input depths z
            method "weightedstairs" means creating a new model that preserves the average values
                use it to resample a depthmodel without loosing its average shape
        :return:
            values of the depthmodel at required depth z
        """
        if interpmethod is None: 
            #no user specification, use the default interpolation method
            interpmethod = self.interpmethod

        if interpmethod == "stairs":
            Z, V = self.stairs()
            return np.interp(xp = Z, fp = V, x = z, left = V[0], right = V[-1])
        elif interpmethod == "weightedstairs":
            # determine the stairs shape of the model
            Z, V = self.stairs()
            Z[-1] = 10. * Z[-2]

            # compute model integral
            U = discrete_time_primitive(Z, V) + V[0]
            u = np.interp(xp = Z, fp = U, x = z, left = U[0], right = U[-1])

            # derivate the integral at required location
            newv = np.concatenate(((u[1:] - u[:-1]) / (z[1:] - z[:-1]),[V[-1]]))

            return newv

        else: raise ValueError('interpmethod %s not valid for %s' % (interpmethod, type(self)))

    #---------------------------------------------
    def copy(self): return copy.deepcopy(self)

    #----------------------------------------------------
    def split(self, thickness):
        """tries to split layers with a resolution of thickness
           while preserving the interfaces already in place
           overwrites self
        """
        assert thickness > 0.
        thck, values = [], []

        for i, (h, v) in enumerate(zip(self.thickness(), self.values)):
            if np.isinf(h): 
                thck.append(np.inf)
                values.append(v)
            else:
                #q, r = diveucl(h, thickness)
                q, r = h // thickness, h % thickness
                if r > 1.e-6:
                    values.append(v)
                    thck.append(r)

                if q > 0.:
                    for j in xrange(int(q)):
                        values.append(v)
                        thck.append(thickness)

        ztop = np.concatenate(([0.], np.cumsum(thck[:-1])))
        ztop   = np.array(ztop)
        values = np.array(values)

        self.__init__(ztop, values)

    #----------------------------------------------------
    def add_layer(self, thickness):
        """add one layer on top of half-space
           overwrites self
        """
        if hasattr(thickness, "__iter__"):
            ztop   = np.concatenate((self.z[:],      self.z[-1] +  np.cumsum(thickness)))
            values = np.concatenate((self.values[:], self.values[-1] * np.ones(len(thickness), float)))
            self.__init__(ztop, values)            
        else:
            ztop   = np.concatenate((self.z[:],      [self.z[-1] + thickness]))
            values = np.concatenate((self.values[:], [self.values[-1]]))
            self.__init__(ztop, values)

    #----------------------------------------------------
    def smooth(self):
        """split every layer in two peace
           overwrites self
        """
        z = np.concatenate((self.z, self.zmid()[:-1], [self.z[-1] * 1.05]))
        v = np.concatenate((self.values, self.values))
        I = np.argsort(z)
        z, v = z[I], v[I]
        #plt.figure()
        #self.show(gca())

        v[0] -= 0.25 * (v[2] - v[1])
        v[2:-1:2] -= 0.5 * (v[2:-1:2] - v[1:-2:2])
        newztop = np.concatenate(([z[0]], 0.5 * (z[1:] + z[:-1]), [z[-1]]))
        newv    = np.interp(newztop, xp = z, fp = v)
        newv[1:-2] = v[1:-1]

        self.__init__(newztop, newv)
        #self.show(gca(), '+-')
        #showme()


        #z = np.unique(np.sort(np.concatenate((self.ztop(), self.zmid()[:-1], self.zbot()[:-1]))))
        #v = self.interp(z)
        #v[0] -= 0.25 * (v[2] - v[1])
        #v[2:-1:2] -= 0.5 * (v[2:-1:2] - v[1:-2:2])
        #newztop = np.concatenate(([z[0]], 0.5 * (z[1:] + z[:-1]), [z[-1]]))
        #newv    = np.interp(newztop, xp = z, fp = v)
        #newv[1:-2] = v[1:-1]
        #print newztop[1:] - newztop[:-1]
        #self.__init__(newztop, newv)       

    # ----------------------------------------------------
    def blur(self, thickness):
        """
        apply gaussian smoothing to the depthmodel
        try to preserve surface and half-space velocity
        overwrites self
        thickness = smoothing thickness in km
        """
        # 1/ interpolate the depthmodel on a linear grid
        pad = 256
        z = np.linspace(0.0, self.z.max(), 1792)
        v = self.interp(z, interpmethod = "stairs")

        # 2/ remove the trend, keep trend coefficients
        a = (v[-1] - v[0]) / (z[-1] - z[0])
        b = v[0]
        v1 = v - (a * z + b)
        v1 = np.concatenate((v1, np.zeros(pad, float)))

        # 3/ perform gaussian blur in spectral domain
        N, dz = len(v1), z[1] - z[0]
        fn = fftfreq(N, 1. / N)

        dn = thickness / dz
        I = np.exp(-0.5 * (((abs(fn) - 0.) / dn) ** 2.)) / (np.sqrt(2. * np.pi) * dn)
        v2 = np.real(ifft(fft(v1) * fft(I)))

        # 4/ taper the edges to preserve surface and half-space amplitude
        v2 = v2[:len(z)] #remove pad
        v2 *= cosTaperwidth(v2, sampling_rate = 1. / dz, width=thickness)

        # 5/ restore trend
        v3 = (a * z + b) + v2

        # #qc
        # plt.figure()
        # gca().plot(z, v, z, v1[:len(z)], z, v2, z, v3)
        # showme()

        # 6/ re-interpolate smoothed version of the model at self's depths
        zi = self.zmid()
        zi[-1] = self.z[-1]
        self.values = np.interp(xp = z, fp = v3, x = zi)

    #----------------------------------------------------
    def __call__(self, z):
        raise Exception('obsolet, use depthmodel1D.interp instead')
        zz, vv = self.stairs()
        return np.interp(xp = zz, fp = vv, x = z)
#----------------------------------------------------
class depthmodel(object):
    """combination of 3 depthmodel1D for vp, vs and density
    the depth grid must be the same for all three models
    """
    def __init__(self, vp, vs, rh, check = True):
        if check:
            assert isinstance(vp, depthmodel1D)
            assert isinstance(vs, depthmodel1D)
            assert isinstance(rh, depthmodel1D)
            assert len(vp) == len(vs) == len(rh)
            assert np.all(vp.z == vs.z)
            assert np.all(vp.z == rh.z)
            assert np.all(vp.values > vs.values)
            assert np.all(vs.values > 0.)
            assert np.all(rh.values > 0.)

        self.vp, self.vs, self.rh = vp, vs, rh

    # ----------------------------------------------------
    def ztop(self): return self.vs.ztop()
    def zbot(self): return self.vs.zbot()
    def zmid(self): return self.vs.zmid()
    def thickness(self): return self.vs.thickness()
    def __len__(self): return len(self.vs)

    #----------------------------------------------------
    def simplify(self):
        #group layers with same values, return it as a new object
        ztop, vp, vs, rh = [self.vp.z[0]], [self.vp.values[0]], [self.vs.values[0]], [self.rh.values[0]]
        for i in xrange(1, len(self.vp)):
            if self.vp.values[i] != self.vp.values[i-1] or \
               self.vs.values[i] != self.vs.values[i-1] or \
               self.rh.values[i] != self.rh.values[i-1]:
                ztop.append(self.vp.z[i])
                vp.append(self.vp.values[i])
                vs.append(self.vs.values[i])
                rh.append(self.rh.values[i])
        return depthmodel(\
            depthmodel1D(ztop, vp),
            depthmodel1D(ztop, vs),
            depthmodel1D(ztop, rh))
                
    #----------------------------------------------------
    def vp_over_vs(self):
        return depthmodel1D(self.ztop(), self.vp.values / self.vs.values)

    #----------------------------------------------------
    def pr(self):return self.vp_over_vs()

    #----------------------------------------------------
    def mu(self):
        return depthmodel1D(self.ztop(), self.rh.values * self.vs.values ** 2.)

    #----------------------------------------------------
    def lamda(self):
        return depthmodel1D(self.ztop(), self.rh.values * (self.vp.values ** 2. - 2. * self.vs.values ** 2.))

    #----------------------------------------------------
    def __str__(self):
        """
        :return: string at mod96 file format
        """
        strout  = 'MODEL.01\nwhatever\nISOTROPIC\nKGS\nFLAT EARTH\n1-D\nCONSTANT VELOCITY\n'
        strout += 'LINE08\nLINE09\nLINE10\nLINE11\n'
        strout += '      H(KM)   VP(KM/S)   VS(KM/S) RHO(GM/CC)     QP         QS       ETAP       ETAS      FREFP      FREFS\n'
        fmt     = "      %.6f    %.6f       %.6f     %.6f           0.         0.       0.         0.        1.         1.\n"
        H = self.thickness()
        H[-1] = 0.

        for h, vp, vs, rh in zip(H, self.vp.values, self.vs.values, self.rh.values):
            strout += fmt % (h, vp, vs, rh)
        return strout

    #----------------------------------------------------
    def write96(self, filename, overwrite = False):
        """
        :param filename: string, name of the ascii file to create at .mod96 format
        :param overwrite: bool, overwrite if exists
        :return: None
        """
        if not overwrite and os.path.exists(filename): raise Exception('%s already exists' % filename)
        with open(filename, 'w') as fid:
            fid.write(self.__str__())

    #---------------------------------------------
    def copy(self): return copy.deepcopy(self)

    #----------------------------------------------------
    def split(self, *args, **kwargs):
        for d in self.vp, self.vs, self.rh:
            d.split(*args, **kwargs)

    #----------------------------------------------------
    def add_layer(self, *args, **kwargs):
        for d in self.vp, self.vs, self.rh:
            d.add_layer(*args, **kwargs)

    #----------------------------------------------------
    def interp(self, ztop, **kwargs):
        vp = depthmodel1D(ztop, self.vp.interp(ztop, **kwargs))
        vs = depthmodel1D(ztop, self.vs.interp(ztop, **kwargs))
        rh = depthmodel1D(ztop, self.rh.interp(ztop, **kwargs))
        return depthmodel(vp, vs, rh)

    #----------------------------------------------------
    def blur(self, thickness):
        for d in self.vp, self.vs, self.rh:
            d.blur(thickness)

    #----------------------------------------------------
    def ztopvpvsrh(self):
        return self.ztop(), self.vp.values, self.vs.values, self.rh.values
    #----------------------------------------------------
    def show(self, ax, *args, **kwargs):

        hdls = []
        kwargs['color'] = "b"
        hdls.append(self.vp.show(ax, *args, **kwargs))
        kwargs['color'] = "g"
        hdls.append(self.vs.show(ax, *args, **kwargs))
        kwargs['color'] = "r"
        hdls.append(self.rh.show(ax, *args, **kwargs))
        return hdls
#----------------------------------------------------
class depthmodel_from_mod96_txt(depthmodel):
    """initiate from a string containing a depthmodel at mod96 format"""
    def __init__(self, mod96txt):
        L = mod96txt.split('\n')
        while True:
            if "H(KM)" in L.pop(0): break
        
        H, VP, VS, RH = [], [], [], []
        while len(L):
            l = L.pop(0)
            l = l.split('/n')[0].split()
            if not len(l): continue
            if len(l) != 10: raise Exception('%s not understood' % l)
            h, vp, vs, rh = np.array(l[:4], float)
            H.append(h); VP.append(vp); VS.append(vs); RH.append(rh)
            
        H, VP, VS, RH = [np.array(x) for x in H, VP, VS, RH]
        ztop = np.concatenate(([0.], H[:-1].cumsum()))

        vp = depthmodel1D(ztop, VP)
        vs = depthmodel1D(ztop.copy(), VS)
        rh = depthmodel1D(ztop.copy(), RH)
        depthmodel.__init__(self, vp, vs, rh)
#----------------------------------------------------
class depthmodel_from_mod96(depthmodel):
    """initiate from a file at mod96 format to read"""
    def __init__(self, filename):
        with open(filename, 'r') as fid:
            while True:
                if "H(KM)" in fid.readline(): break
            
            H, VP, VS, RH = [], [], [], []
            while True:
                l = fid.readline()
                if l == "": break
                l = l.split('/n')[0].split()
                if not len(l): continue
                if len(l) != 10: raise Exception('%s not understood' % l)
                h, vp, vs, rh = np.array(l[:4], float)
                H.append(h); VP.append(vp); VS.append(vs); RH.append(rh)
            
        H, VP, VS, RH = [np.array(x) for x in H, VP, VS, RH]
        ztop = np.concatenate(([0.], H[:-1].cumsum()))

        vp = depthmodel1D(ztop, VP)
        vs = depthmodel1D(ztop.copy(), VS)
        rh = depthmodel1D(ztop.copy(), RH)
        depthmodel.__init__(self, vp, vs, rh)
#----------------------------------------------------
class depthpdf(object):
    def __init__(self, z, v):
        assert np.all(z[1:] > z[:-1])
        assert np.all(v[1:] > v[:-1])
        self.z = z
        self.v = v
        self.H = np.zeros((len(z), len(v)), float)

    #------------------------------------------------
    def write(self, filename):
        assert filename.split('.')[-1] == "zpdf"
        pkl((self.z, self.v, self.H), filename)

    #------------------------------------------------
    def append(self, m1D, **kwargs):
        #kwargs are passed to m1D.interp 
        self.appendN(m1D, Ntimes = 1, **kwargs)

    #------------------------------------------------
    def appendN(self, m1D, Ntimes=1, **kwargs):
        """append the same model Ntimes times in the histogram"""
        def search_sorted_nearest(a, v, lena):
            n = np.searchsorted(a, v)
            if n == lena: return lena - 1
            elif n == 0: return 0
            else:
                vsup = a[n]
                vinf = a[n-1]
                if abs(vsup - v) <= abs(v - vinf):
                    return n
                return n-1


        zstairs, vstairs = m1D.stairs()
        zstairs[-1] = self.z.max() + 1.#np.max([self.z.max(), 1.5 * zstairs[-2] + 1.0])
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
            zmin, zmax = zstairs[i], zstairs[i+1]
            i0 = search_sorted_nearest(self.z, zmin, lenz)#np.clip(np.searchsorted(self.z, zmin), 0, len(self.z) - 1)
            i1 = search_sorted_nearest(self.z, zmax, lenz) + 1#np.clip(np.searchsorted(self.z, zmax), 0, len(self.z))
            jj = search_sorted_nearest(self.v, vstairs[i], lenv) #np.clip(np.searchsorted(self.v, vstairs[i]), 0,len(self.v) - 1)
            self.H[i0:i1,jj] += float(Ntimes)
    

    #------------------------------------------------
    def show(self, ax, **kwargs):
        I = np.zeros(len(self.z))
        for i in xrange(len(self.z)):
            I[i] = discrete_time_primitive(self.v, self.H[i, :], area = True)
        I = discrete_time_primitive(self.z, I, area = True)


        H = np.ma.masked_where(self.H == 0., self.H)
        v_ = np.concatenate((self.v, [self.v[-1] + self.v[-1] - self.v[-2]]))
        z_ = np.concatenate((self.z, [self.z[-1] + self.z[-1] - self.z[-2]]))
        ax.pcolormesh1(v_, z_, H / I, **kwargs)
        if not ax.yaxis_inverted():
            ax.invert_yaxis()
    #------------------------------------------------
    def purcentile(self, p):
        assert 0. <= p <= 1.
        P = np.zeros_like(self.z) * np.nan
        for i, z in enumerate(self.z):
            #x = self.v
            y = np.cumsum(self.H[i, :])
            ymin, ymax = y.min(), y.max()
            if ymax > ymin:
                y = (y - y.min()) / (y.max() - y.min())
                P[i] = np.interp(p, xp=y, fp=self.v)
            elif ymax == ymin:
                P[i] = self.v[0]
        assert not np.isnan(P).any()
        return self.z, P

#------------------------------------------------
class depthpdf_from_zpdffile(depthpdf):
    def __init__(self, filename):
        self.z, self.v, self.H = unpkl(filename)
#------------------------------------------------
def dmstats(dms, percentiles = [0.16, 0.5, 0.84], Ndepth = 100, Nvalue = 100, weights = None):
    assert np.all([isinstance(dm, depthmodel) for dm in dms])
    assert np.all([0 < p < 1 for p in percentiles])
    assert len(percentiles) == len(np.unique(percentiles))
    if weights is None: weights = np.ones(len(dms))
    else: assert len(weights) == len(dms)

    zmax = -np.inf
    vsmin, vsmax = np.inf, -np.inf
    vpmin, vpmax = np.inf, -np.inf
    rhmin, rhmax = np.inf, -np.inf
    prmin, prmax = np.inf, -np.inf
    for dm in dms:
        zmax  = np.max([zmax, 1.1 * dm.vs.z[-1]])
        vsmin = np.min([vsmin, dm.vs.values.min()])
        vsmax = np.max([vsmax, dm.vs.values.max()])
        vpmin = np.min([vpmin, dm.vp.values.min()])
        vpmax = np.max([vpmax, dm.vp.values.max()])
        rhmin = np.min([rhmin, dm.rh.values.min()])
        rhmax = np.max([rhmax, dm.rh.values.max()])
        prmin = np.min([prmin, dm.pr().values.min()])
        prmax = np.max([prmax, dm.pr().values.max()])

    zbins = np.linspace(0., zmax, Ndepth)
    vspdf = depthpdf(z = zbins, v = np.linspace(vsmin, vsmax, Nvalue)) 
    vppdf = depthpdf(z = zbins, v = np.linspace(vpmin, vpmax, Nvalue))
    rhpdf = depthpdf(z = zbins, v = np.linspace(rhmin, rhmax, Nvalue))
    prpdf = depthpdf(z = zbins, v = np.linspace(prmin, prmax, Nvalue))


    for dm, weight in zip(dms, weights):
        vspdf.appendN(dm.vs, Ntimes=weight)
        vppdf.appendN(dm.vp, Ntimes=weight)
        rhpdf.appendN(dm.rh, Ntimes=weight)
        prpdf.appendN(dm.pr(), Ntimes=weight)

    for p in percentiles:
        zpc, vspc = vspdf.purcentile(p)
        _,   vppc = vppdf.purcentile(p)
        _,   rhpc = rhpdf.purcentile(p)
        _,   prpc = prpdf.purcentile(p)
        vppc = depthmodel1D(zpc, vppc)#.simplify()
        vspc = depthmodel1D(zpc, vspc)#.simplify()
        rhpc = depthmodel1D(zpc, rhpc)#.simplify()
        prpc = depthmodel1D(zpc, prpc)#.simplify()

        #dmpc = depthmodel(\).simplify()
        #yield p, dmpc #, (zmed, vpmed, vsmed, rhmed, prmed)
        yield p, (vppc, vspc, rhpc, prpc)

# ------------------------------------------------
class UserStacker(object):
    def __init__(self, zmax, Ndepth, vpmin, vpmax, vsmin, vsmax, rhmin, rhmax, prmin, prmax, Nvalue):
        zbins = np.linspace(0., zmax, Ndepth)
        self.vspdf = depthpdf(z=zbins, v=np.linspace(vsmin, vsmax, Nvalue))
        self.vppdf = depthpdf(z=zbins, v=np.linspace(vpmin, vpmax, Nvalue))
        self.rhpdf = depthpdf(z=zbins, v=np.linspace(rhmin, rhmax, Nvalue))
        self.prpdf = depthpdf(z=zbins, v=np.linspace(prmin, prmax, Nvalue))

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
# ------------------------------------------------
def dmstats1(dms, percentiles=[0.16, 0.5, 0.84], Ndepth=100, Nvalue=100, weights=None, **mapkwargs):
    assert np.all([isinstance(dm, depthmodel) for dm in dms])
    assert np.all([0 < p < 1 for p in percentiles])
    assert len(percentiles) == len(np.unique(percentiles))
    if weights is None: weights = np.ones(len(dms))
    else: assert len(weights) == len(dms)

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
    s0 = UserStacker(zmax, Ndepth, vpmin, vpmax, vsmin, vsmax, rhmin, rhmax, prmin, prmax, Nvalue) #the stacker to be reproduced in each independent workspace (deep copy)
    with StackAsync(s0, JobGen(), Verbose = False, **mapkwargs) as sa:
        S = UserStacker(zmax, Ndepth, vpmin, vpmax, vsmin, vsmax, rhmin, rhmax, prmin, prmax, Nvalue) #initiate the final stack
        for jobids, (s, wname), Tgen, Tpro in sa: #receive partial stacks
            sa.communicate("Stacker %s stacked %6d jobs in %.6fs" % (wname, len(jobids), Tgen + Tpro))
            S += s #merge all partial stacks using method __iadd__

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


# class depthmodel_from_vs_pr_rh(depthmodel):
#     """
#     same as depthmodel but
#     """
#     def __init__(self, vs, pr, rh, check = True):
#         if check:
#             assert isinstance(vs, depthmodel1D)
#             assert isinstance(pr, depthmodel1D)
#             assert isinstance(rh, depthmodel1D)
#             assert len(vs) == len(pr) == len(rh)
#             assert np.all(vs.z == pr.z)
#             assert np.all(vs.z == rh.z)
#             assert np.all(vs.values > 0.)
#             assert np.all(pr.values > 1.)
#             assert np.all(rh.values > 0.)
#
#         self.vs, self.pr, self.rh = vs, pr, rh
#     #----------------------------------------------------
#     def vp(self):
#         return depthmodel1D(self.ztop(), self.pr.values * self.vs.values)
#     #----------------------------------------------------
#     def vp_over_vs(self): return self.pr
#     #----------------------------------------------------
#     def lamda(self):
#         return depthmodel1D(self.ztop(), self.rh.values * (self._vp.values ** 2. - 2. * self.vs.values ** 2.))
#
#     #----------------------------------------------------
#     def __str__(self):
#         strout  = 'MODEL.01\nwhatever\nISOTROPIC\nKGS\nFLAT EARTH\n1-D\nCONSTANT VELOCITY\n'
#         strout += 'LINE08\nLINE09\nLINE10\nLINE11\n'
#         strout += '      H(KM)   VP(KM/S)   VS(KM/S) RHO(GM/CC)     QP         QS       ETAP       ETAS      FREFP      FREFS\n'
#         fmt     = "      %.6f    %.6f       %.6f     %.6f           0.         0.       0.         0.        1.         1.\n"
#         H = self.thickness()
#         H[-1] = 0.
#
#         for h, vp, vs, rh in zip(H, self.vp().values, self.vs.values, self.rh.values):
#             strout += fmt % (h, vp, vs, rh)
#         return strout
#
#     #----------------------------------------------------
#     def split(self, *args, **kwargs):
#         for d in self.vs, self.pr, self.rh:
#             d.split(*args, **kwargs)
#
#     #----------------------------------------------------
#     def add_layer(self, *args, **kwargs):
#         for d in self.vs, self.pr, self.rh:
#             d.add_layer(*args, **kwargs)
#
#     #----------------------------------------------------
#     def interp(self, ztop, **kwargs):
#         vs = depthmodel1D(ztop, self.vs.interp(ztop, **kwargs))
#         pr = depthmodel1D(ztop, self.pr.interp(ztop, **kwargs))
#         rh = depthmodel1D(ztop, self.rh.interp(ztop, **kwargs))
#         return depthmodel_from_vs_pr_rh(vs, pr, rh)
#
#     #----------------------------------------------------
#     def show(self, ax, *args, **kwargs):
#         kwargs['color'] = "g"
#         self.vs.show(ax, *args, **kwargs)
#         kwargs['color'] = "k"
#         self.pr.show(ax, *args, **kwargs)
#         kwargs['color'] = "r"
#         self.rh.show(ax, *args, **kwargs)
#####################################################
# def bezier2D(xcontrol, ycontrol, Nbezier):
#     t = np.linspace(0., 1., Nbezier)
#     N = len(xcontrol)
#     bernsteins = np.array([comb(N - 1, i) for i in xrange(N)])
#     xc, yc = xcontrol * bernsteins, ycontrol * bernsteins
#     x, y = np.zeros(len(t), float), np.zeros(len(t), float)
#     for i in np.arange(N):
#         x += xc[i] * (t ** i) * ((1. - t) ** (N - 1. - i))
#         y += yc[i] * (t ** i) * ((1. - t) ** (N - 1. - i))
#     return x, y
# #----------------------------------------------------
# def bezier4points(xa, ya, xb, yb, ra, ta, rb, tb):
#     x0, y0 = xa, ya
#     x1, y1 = xa + ra * np.sin(ta * np.pi / 180.), ya + ra * np.cos(ta * np.pi / 180.)
#     x2, y2 = xb + rb * np.sin(tb * np.pi / 180.), yb + rb * np.cos(tb * np.pi / 180.)
#     x3, y3 = xb, yb
#     X, Y = bezier2D([x0, x1, x2, x3], [y0, y1, y2, y3], Nbezier = 100)
#     return X, Y, [x0, x1, x2, x3], [y0, y1, y2, y3]
# #----------------------------------------------------
# class cubicbezier(object):
#     def __init__(self, zcontrol, vcontrol, rcontrol = None, tcontrol = None):
#         """r = longueur du segment de pente
#            t = angle mesure au point de control, compte positivement dans le sens anti horaire depuis la verticale"""
#
#         self.zcontrol, self.vcontrol = [np.array(v) for v in zcontrol, vcontrol]
#         assert np.all(self.zcontrol[1:] > self.zcontrol[:-1])
#
#         #---------------
#         if tcontrol is None:
#             self.tcontrol = np.zeros_like(self.vcontrol)
#             self.tcontrol[1:-1] = np.arctan((self.vcontrol[2:] - self.vcontrol[:-2]) / (self.zcontrol[2:] - self.zcontrol[:-2])) * 180. / np.pi
#         else:
#             self.tcontrol = np.array(tcontrol)
#             assert len(self.tcontrol) == len(self.vcontrol)
#         assert np.all(self.tcontrol >= -90.) and np.all(self.tcontrol <= 90.)
#         #---------------
#         if rcontrol is None:
#             self.rcontrol = np.nan * np.ones_like(self.vcontrol)
#             dz = self.zcontrol[1:] - self.zcontrol[:-1]
#             self.rcontrol[1:-1] = np.min(np.concatenate(([dz[1:]], [dz[:-1]]), axis = 0), axis = 0) / 2.
#             self.rcontrol[0] = dz[0] / 2.
#             self.rcontrol[-1] = dz[-1] / 2.
#         else:
#             self.rcontrol = np.array(rcontrol)
#             assert len(self.rcontrol) == len(self.vcontrol)
#
#         #---------------
#         self.Z, self.V = [], []
#         for j in xrange(len(self.zcontrol) - 1):
#             vj, zj = bezier4points(\
#                 xa = self.vcontrol[j],       ya = self.zcontrol[j],
#                 xb = self.vcontrol[j + 1],   yb = self.zcontrol[j + 1],
#                 ra = self.rcontrol[j],       rb = self.rcontrol[j + 1],
#                 ta = self.tcontrol[j],       tb = self.tcontrol[j + 1] + 180.)[:2]
#             self.Z += list(zj[:-1])
#             self.V += list(vj[:-1])
#         self.Z.append(zj[-1])
#         self.V.append(vj[-1])
#         self.Z, self.V = np.array(self.Z).flat[:], np.array(self.V).flat[:]
#         #if np.any(self.Z[1:] < self.Z[:-1]): raise ParametrizationError('')
#
#     #------------------------------------------------
#     def plot(self, ax, **kwargs):
#         ax.plot(self.V, self.Z, **kwargs)
#
#     #------------------------------------------------
#     def plotcontrols(self, ax, **kwargs):
#         for j in xrange(len(self.zcontrol) - 1):
#             vj, zj, (x0, x1, x2, x3), (y0, y1, y2, y3) = bezier4points(\
#                 xa = self.vcontrol[j],       ya = self.zcontrol[j],
#                 xb = self.vcontrol[j + 1],   yb = self.zcontrol[j + 1],
#                 ra = self.rcontrol[j],       rb = self.rcontrol[j + 1],
#                 ta = self.tcontrol[j],       tb = self.tcontrol[j + 1] + 180.)
#             ax.plot([x0, x1, x2, x3], [y0, y1, y2, y3], 'o', **kwargs)
#             ax.plot([x0, x1], [y0, y1], '--', **kwargs)
#             ax.plot([x2, x3], [y2, y3], '--', **kwargs)
# #----------------------------------------------------
# class depthmodel1D_from_cubicbezier(depthmodel1D):
#     def __init__(self, ztop, zbezier, vbezier, rbezier=None, tbezier=None, vmin = -np.inf, vmax = +np.inf, interpmethod = "bezier"):
#         assert interpmethod in ['stairs', 'bezier']
#         c = cubicbezier(zbezier, vbezier, rbezier, tbezier)
#         z, v = c.Z, c.V
#
#         assert np.all((vmin < v) & (v < vmax))
#         assert np.all((z[1:] - z[:-1]) > 0.)
#         assert ztop[0] == zbezier[0] == 0.
#         assert abs(ztop[-1] - zbezier[-1]) < ztop[-1] / 10000.
#
#
#         method = 1
#         if method == 1:
#             """
#             convert bezier curve into a stairs model
#             mid layers = bezier curve
#             """
#             ztop = np.array(ztop, float)
#             zmid = np.concatenate((ztop[:-1] + 0.5 * (ztop[1:] - ztop[:-1]), [1.01 * ztop[-1]]))
#             vmid = np.interp(x = zmid, xp = z, fp = v)
#             depthmodel1D.__init__(self, ztop = ztop, values = vmid)
#         elif method == 2:
#             """
#             convert bezier curve into a stairs model
#             equal area, make sure the Nbezier parameter is high enough
#             """
#             vv    = discrete_time_primitive(z, v, area = False)
#             vvbot = np.interp(x = ztop[1:],  xp = z, fp = vv)
#             vvtop = np.interp(x = ztop[:-1], xp = z, fp = vv)
#             vvv = (vvbot - vvtop) / (ztop[1:] - ztop[:-1])
#             vvv = np.concatenate((vvv, [v[-1]]))
#             depthmodel1D.__init__(self, ztop = ztop, values = vvv)
#
#         self._cubicbezier = c
#         self.interpmethod = interpmethod
#
#     def interp(self, z, interpmethod = None):
#         if interpmethod is None: interpmethod = self.interpmethod
#         if interpmethod == "stairs":
#             return depthmodel1D.interp(self, z, interpmethod = "stairs")
#         elif interpmethod == "bezier":
#             return np.interp(xp = self._cubicbezier.Z, fp = self._cubicbezier.V, x = z, left = self._cubicbezier.V[0], right = self._cubicbezier.V[-1])
#         else: raise ValueError('interpmethod %s not valid for %s' % (interpmethod, type(self)))
#
#     def show(self, ax, *args, **kwargs):
#         depthmodel1D.show(self, ax, *args, **kwargs)
#         self._cubicbezier.plot(ax, color = [.5, .5, .5])
#         #self._cubicbezier.plotcontrols(ax, color = [.5, .5, .5])
# #----------------------------------------------------
# def example_depthmodel1D_from_cubicbezier():
#
#
#     ztop = depthspace(3., 7)
#     d = depthmodel1D_from_cubicbezier(\
#         zbezier = [0.,  1.,  3.],
#         vbezier = [2.0, 1.,  2.2],
#         rbezier = [0.5, 0.5, 0.5],
#         tbezier = [90., 45., 90.],
#         ztop    = ztop,
#         method  = 1)
#     d._cubicbezier.plot(gca(), color = 'r')
#     d._cubicbezier.plotcontrols(gca(), color = 'r')
#     d.show(gca(), 'b.-')
#
#
#     d = depthmodel1D_from_cubicbezier(\
#         zbezier = [0.,  1.,  3.],
#         vbezier = [2.0, 1.,  2.2],
#         rbezier = [0.5, 0.5, 0.5],
#         tbezier = [90., 45., 90.],
#         ztop    = ztop,
#         method  = 2)
#     d._cubicbezier.plot(gca(), color = 'r')
#     d._cubicbezier.plotcontrols(gca(), color = 'r')
#     d.show(gca(), 'gx-')
#
#     gca().set_aspect(1.)
#     showme()
#####################################################