from srfpython.depthdisp.mod96 import packmod96, unpackmod96
from srfpython.utils import discrete_time_primitive, cosTaperwidth
from scipy.fftpack import fft, ifft, fftfreq
import numpy as np
import copy
import os


# ----------------------------------------------------
def gardner74(vp):
    """Gardner's density law for sedimentary rocks, 1974
       vp between 1.5 and 6.1 km/s
    """
    return 1.74 * vp ** 0.25 #density in g/cm3


# ----------------------------------------------------
def depthspace(zbot, n):
    """
    zbot  = top of the half space
    n     = desired number of layers
    """
    z = np.logspace(np.log10(0.1), np.log10(3.), n)
    z = zbot * (z - z.min()) / (z.max() - z.min())
    return z

# -------------------------------------------------
class depthmodel1D(object):
    """self.z = ztop !!!
       can handle irregular layers
    """

    # -------------------------------------------------
    def __init__(self, ztop, values, interpmethod="stairs"):
        #last value of v stand for half space
        assert len(ztop) == len(values)
        assert np.all(ztop[1:] > ztop[:-1])
        assert ztop[0] == 0.
        self.values, self.z = np.array(values), np.array(ztop)
        self.interpmethod = interpmethod #to be used by default

    # -------------------------------------------------
    def __len__(self):
        return len(self.z)#number of layers including the lower half space

    # -------------------------------------------------
    def ztop(self):
        return self.z

    # -------------------------------------------------
    def thickness(self):
        return np.concatenate((self.z[1:] - self.z[:-1], [np.inf]))

    # -------------------------------------------------
    def zbot(self):
        return self.z + self.thickness()

    # -------------------------------------------------
    def zmid(self):
        return self.z + 0.5 * self.thickness()

    # -------------------------------------------------
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

    # -------------------------------------------------
    def dontshow(self):
        """compute the coordinates of the line for display as in self.show but do not diplay anything"""
        zz, vv = self.stairs()
        zz[-1] = np.max([1.5 * self.z[-1], 3.0])  # cannot plot np.inf...
        return vv, zz

    # -------------------------------------------------
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

    # -------------------------------------------------
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

        else:
            raise ValueError('interpmethod %s not valid for %s' % (interpmethod, type(self)))

    # -------------------------------------------------
    def copy(self):
        return copy.deepcopy(self)

    # -------------------------------------------------
    def split(self, thickness):
        """tries to split layers with a resolution of thickness
           while preserving the interfaces already in place
           overwrites self
        """
        assert thickness > 0.
        # thck, values = [], []
        #
        # for i, (h, v) in enumerate(zip(self.thickness(), self.values)):
        #     if np.isinf(h):
        #         thck.append(np.inf)
        #         values.append(v)
        #     else:
        #         #q, r = diveucl(h, thickness)
        #         q, r = h // thickness, h % thickness
        #         print q, r
        #         if r > 1.e-6:
        #             values.append(v)
        #             thck.append(r)
        #
        #         if q > 0.:
        #             for j in xrange(int(q)):
        #                 values.append(v)
        #                 thck.append(thickness)
        #
        # ztop = np.concatenate(([0.], np.cumsum(thck[:-1])))
        # ztop   = np.array(ztop)
        # values = np.array(values)
        newztop = np.concatenate((self.z, np.arange(0., self.z[-1], thickness)))
        newztop = np.round(newztop, 6)
        newztop = np.sort(np.unique(newztop))

        newvalues = self.interp(newztop+1.e-12, interpmethod="stairs")
        self.__init__(newztop, newvalues)

    # -------------------------------------------------
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

    # -------------------------------------------------
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

    # -------------------------------------------------
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

    # -------------------------------------------------
    def __call__(self, z):
        raise Exception('obsolet, use depthmodel1D.interp instead')
        zz, vv = self.stairs()
        return np.interp(xp = zz, fp = vv, x = z)


# -------------------------------------------------
class depthmodel(object):
    def __init__(self, vp, vs, rh):
        """initiate with 3 depthmodel1D objects, must have the same ztop values"""
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

    # -------------------------------------------------
    def ztop(self):
        return self.vs.ztop()

    # -------------------------------------------------
    def zbot(self):
        return self.vs.zbot()

    # -------------------------------------------------
    def zmid(self):
        return self.vs.zmid()

    # -------------------------------------------------
    def thickness(self):
        return self.vs.thickness()

    # -------------------------------------------------
    def __len__(self):
        return len(self.vs)

    # -------------------------------------------------
    def simplify(self):
        # group layers with same values, return it as a new object
        ztop, vp, vs, rh = [self.vp.z[0]], [self.vp.values[0]], [self.vs.values[0]], [self.rh.values[0]]
        for i in xrange(1, len(self.vp)):
            if self.vp.values[i] != self.vp.values[i - 1] or \
                            self.vs.values[i] != self.vs.values[i - 1] or \
                            self.rh.values[i] != self.rh.values[i - 1]:
                ztop.append(self.vp.z[i])
                vp.append(self.vp.values[i])
                vs.append(self.vs.values[i])
                rh.append(self.rh.values[i])
        return depthmodel( \
            depthmodel1D(ztop, vp),
            depthmodel1D(ztop, vs),
            depthmodel1D(ztop, rh))

    # -------------------------------------------------
    def vp_over_vs(self):
        return depthmodel1D(self.ztop(), self.vp.values / self.vs.values)

    # -------------------------------------------------
    def pr(self):
        return self.vp_over_vs()

    # -------------------------------------------------
    def mu(self):
        return depthmodel1D(self.ztop(), self.rh.values * self.vs.values ** 2.)

    # -------------------------------------------------
    def lamda(self):
        return depthmodel1D(self.ztop(), self.rh.values * (self.vp.values ** 2. - 2. * self.vs.values ** 2.))

    # -------------------------------------------------
    def __str__(self):
        return packmod96(self.vp.z, self.vp.values, self.vs.values, self.rh.values)

    # -------------------------------------------------
    def write96(self, filename):
        with open(filename, 'w') as fid:
            fid.write(self.__str__())

    # -------------------------------------------------
    def copy(self):
        return copy.deepcopy(self)

    # -------------------------------------------------
    def split(self, *args, **kwargs):
        for d in self.vp, self.vs, self.rh:
            d.split(*args, **kwargs)

    # -------------------------------------------------
    def add_layer(self, *args, **kwargs):
        for d in self.vp, self.vs, self.rh:
            d.add_layer(*args, **kwargs)

    # -------------------------------------------------
    def interp(self, ztop, **kwargs):
        vp = depthmodel1D(ztop, self.vp.interp(ztop, **kwargs))
        vs = depthmodel1D(ztop, self.vs.interp(ztop, **kwargs))
        rh = depthmodel1D(ztop, self.rh.interp(ztop, **kwargs))
        return depthmodel(vp, vs, rh)

    # -------------------------------------------------
    def blur(self, thickness):
        for d in self.vp, self.vs, self.rh:
            d.blur(thickness)

    # -------------------------------------------------
    def ztopvpvsrh(self):
        return self.ztop(), self.vp.values, self.vs.values, self.rh.values

    # -------------------------------------------------
    def show(self, ax, *args, **kwargs):

        hdls = []
        kwargs['color'] = "b"
        hdls.append(self.vp.show(ax, label="$ V_p (km/s) $", *args, **kwargs))
        kwargs['color'] = "g"
        hdls.append(self.vs.show(ax, label="$ V_s (km/s) $", *args, **kwargs))
        kwargs['color'] = "r"
        hdls.append(self.rh.show(ax, label=r"$ \rho (g/cm^3) $", *args, **kwargs))
        ax.set_ylabel('$ depth (km) $')
        return hdls

# -------------------------------------------------
class depthmodel_from_arrays(depthmodel):
    def __init__(self, z, vp, vs, rh):
        """initiate with arrays, skip verification for same depth array"""
        assert len(z) == len(vp) == len(vs) == len(rh)
        z, vp, vs, rh = [np.asarray(_, float) for _ in z, vp, vs, rh]
        assert z[0] == 0
        assert np.all(z[1:] > z[:-1])
        assert np.all(vs > 0.)
        assert np.all(rh > 0.)
        assert np.all(vp / vs >= np.sqrt(4 / 3.))

        self.vp, self.vs, self.rh = [depthmodel1D(z, _) for _ in vp, vs, rh]


# -------------------------------------------------
class depthmodel_from_mod96string(depthmodel):
    def __init__(self, mod96string):
        _, Z, _, VP, VS, RHO, _, _, _, _, _, _ = \
            unpackmod96(mod96string)

        depthmodel.__init__(self,
                            vp=depthmodel1D(Z, VP),
                            vs=depthmodel1D(Z, VS),
                            rh=depthmodel1D(Z, RHO))


# -------------------------------------------------
class depthmodel_from_mod96(depthmodel_from_mod96string):
    def __init__(self, filename):
        with open(filename, 'r') as fid:
            depthmodel_from_mod96string.__init__(self, "".join(fid.readlines()))
