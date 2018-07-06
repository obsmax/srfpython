from srfpython.depthdisp.dispcurves import Claw, freqspace
from scipy.fftpack import fftfreq, fft, ifft
from signalfuncs import detrend, taperwidth, gaussbandpass, bandpass
import numpy as np
"""
compute synthetic seismograms or correlation functions
in both cases, the Greens functions determines the time array according to the source and receiver(s) positions
"""


def isfloat(x):
    return isinstance(x, float) and not np.isnan(x) and not np.isinf(x)


# ---------------------------
class Green(object):
    """
    simplistic greens function in a 1D elastic dispersive model
    **  account only for the phase dispersion and the geometrical spreading of
        a multi modal surface wave

    """
    # ---------------------------
    def __init__(self, claws, fmin=0.2, fmax=5., order=4.0, dt=1. / 12.):  # , proj = mp):
        """
        :param claws: list of Claw objects, phase velocity dispersion laws, first one must have interpolationmode = to 0
        :param fmin: lower frequency for butterworth or central freq for gaussian bandpass, (Hz)
        :param fmax: upper frequency for butterworth or central freq for gaussian bandpass, (Hz)
        :param order: filter order for butterworth or gaussian bandpass
        :param dt: sampling interval (s)
        """
        self.claws = claws  #
        self.fmin = fmin  # 0.2
        self.fmax = fmax  # 5.
        self.order = order  # 4.
        self.dt = dt  # 1. / 12.

        assert self.fmax <= 0.5 / self.dt  # nyquist criterion
        assert 0. < self.fmin < self.fmax

        # find lowest group velocities
        self.umin, self.umax = np.inf, -np.inf
        for claw in self.claws:
            ftmp = freqspace(self.fmin, self.fmax, 100, 'flog')
            utmp = claw.to_ulaw()(ftmp)
            utmp = utmp[~np.isnan(utmp)]
            if len(utmp):
                self.umin = np.min([self.umin, utmp.min()])
                self.umax = np.max([self.umax, utmp.max()])
        assert isfloat(self.umin)
        assert isfloat(self.umax)
        assert 0. <= self.umin <= self.umax

        print "group velocity range : %f - %fkm/s" % (self.umin, self.umax)

    # ---------------------------
    def __call__(self, ts, xs, ys, xr, yr, delta=0.05):
        """
        compute synthetic waveform
        :param ts: source time, s
        :param xs: source x coordinate, km
        :param ys: source y coordinate, km
        :param xr: receiver x coordinate, km
        :param yr: receiver y coordinate, km
        :param delta: coefficient for lag time adjustment
        :return starttime: starttime of the trace in s
        :return npts: number of samples in the trace
        :return dt: sampling interval (s)
        :return y: data array
        """
        assert np.all([not hasattr(w, "__iter__") for w in ts, xs, ys, xr, yr])
        # ---------------------------- adjust time array to avoid wrapping effects
        distance = ((xs - xr) ** 2. + (ys - yr) ** 2.) ** 0.5
        tmin = distance / (self.umax * (1. + delta)) - 3. / self.fmin  # first arrival
        tmax = distance / (self.umin / (1. + delta)) + 3. / self.fmin
        if True:  # align sampling to grid
            tmin = np.floor((ts + tmin) / self.dt) * self.dt - ts

        npts = int(np.ceil((tmax - tmin) / self.dt))
        if npts >= 360000:
            print tmin, tmax, npts
            raise Exception('npts > 360000')

        # ----------------------------
        starttime = ts + tmin

        # ----------------------------
        nu = fftfreq(npts, self.dt)
        Y = np.zeros(npts, complex)
        for claw in self.claws:  # sum over modal contributions
            c = claw(nu)
            I = ~np.isnan(c)
            if I.any():
                Y[I] += np.exp(-2.j * np.pi * nu[I] * (distance / c[I] - tmin))
            else:
                print "warning", claw
        # ----------------------------
        if distance: Y /= np.sqrt(distance)  # geometrical spreading

        y = np.real(ifft(Y))

        y = detrend(y)
        y *= taperwidth(y, 1. / self.dt, .5 / self.fmin)
        if self.fmin == self.fmax:
            y = gaussbandpass(y, 1. / self.dt, self.fmin, self.order)
        else:
            y = bandpass(y, 1. / self.dt, self.fmin, self.fmax, self.order, True)
        return starttime, npts, self.dt, y

    # ---------------------------
    def ccf(self, xS, yS, xA, yA, xB, yB,
            delta=0.05,
            derivate=True,
            geometrical_spreading=True):
        """compute synthetic cross correlation functions for a pointwise source
        adjusts the lag time automatically to prevent wrapping

        :param xS: source x coordinate, km
        :param yS: source y coordinate, km
        :param xA: A station x coordinate, km
        :param yA: A station y coordinate, km
        :param xB: B station x coordinate, km
        :param yB: B station x coordinate, km
        :param delta: coefficient for lag time adjustment
        :param derivate: if you want to derivate the correlation functions
        :param geometrical_spreading: if you want to include geometrical attenuations
        :return starttime: starttime of the trace in s
        :return npts: number of samples in the trace
        :return dt: sampling interval (s)
        :return y: data array
        """

        assert np.all([not hasattr(w, "__iter__") for w in xS, yS, xA, yA, xB, yB])
        # if switch(lon1=xA, lat1=yA, lon2=xB, lat2=yB, chnl1=None, chnl2=None):
        #     return self.ccf(xS, yS,
        #                     xA=xB, yA=yB, zA=zB,
        #                     xB=xA, yB=yA, zB=zA,
        #                     networkA=networkB,
        #                     stationA=stationB,
        #                     locationA=locationB,
        #                     networkB=networkA,
        #                     stationB=stationA,
        #                     locationB=locationA,
        #                     derivate=derivate, **kwargs)
        # ----------------------
        SA = ((xS - xA) ** 2. + (yS - yA) ** 2.) ** 0.5
        SB = ((xS - xB) ** 2. + (yS - yB) ** 2.) ** 0.5
        distance = ((xA - xB) ** 2. + (yA - yB) ** 2.) ** 0.5

        if SA >= SB:
            #            tmin = (SA - SB) / self.umax - 3. / self.fmin
            #            tmax = (SA - SB) / self.umin + 3. / self.fmin
            tmin = (SA - SB) / (self.umax * (1. + delta)) - 3. / self.fmin
            tmax = (SA - SB) / (self.umin / (1. + delta)) + 3. / self.fmin
        else:
            tmin = (SA - SB) / (self.umin / (1. + delta)) - 3. / self.fmin
            tmax = (SA - SB) / (self.umax * (1. + delta)) + 3. / self.fmin

        tmin = np.floor(tmin / self.dt) * self.dt
        npts = int(np.ceil((tmax - tmin) / self.dt))
        if npts % 2: npts += 1
        tmax = tmin + npts * self.dt

        # ----------------------------
        starttime = tmin

        # ----------------------------
        nu = fftfreq(npts, self.dt)
        Y = np.zeros(npts, complex)
        for nclaw in self.claws:  # sum over modal contributions
            cn = nclaw(nu)
            In = ~np.isnan(cn)
            for mclaw in self.claws:  # sum over modal contributions
                cm = nclaw(nu)
                Im = ~np.isnan(cm)
                I = In & Im
                if I.any():
                    # computes uA * conj(uB)
                    Y[I] += np.exp(-2.j * np.pi * nu[I] * (SA / cn[I] - SB / cm[I] - tmin))

        if derivate: Y *= 2.j * np.pi * nu
        if SA and geometrical_spreading: Y /= float(np.sqrt(SA))
        if SB and geometrical_spreading: Y /= float(np.sqrt(SB))
        y = np.real(ifft(Y))
        # ----------------------------
        y = detrend(y)
        y *= taperwidth(y, 1. / self.dt, .5 / self.fmin)
        if self.fmin == self.fmax:
            y = gaussbandpass(y, 1. / self.dt, self.fmin, self.order)
        else:
            y = bandpass(y, 1. / self.dt, self.fmin, self.fmax, self.order, True)

        return starttime, npts, self.dt, y


# ##################################################
if __name__ == "__main__":

    from srfpython.standalone.display import *
    from testlaws import laws

    if True:
        plt.figure()
        G = Green(laws, fmin=0.2, fmax=2., order=4.0, dt=1. / 12.)

        ts, xs, ys = 0., 0., 0.
        for xr in np.arange(1., 20., .25):
            yr = 0.
            starttime, npts, dt, y = G(ts, xs, ys, xr, yr)
            d = np.sqrt((xs - xr) ** 2. + (ys - yr) ** 2.)
            t = starttime + np.arange(npts) * dt
            gca().plot(t, 4. * y + d, "k", alpha=.3)
        gca().set_xlabel('time (s)')
        gca().set_ylabel('distance (km)')
        showme()

    if True:
        plt.figure(figsize=(12, 6))
        G = Green(laws, dt=1. / 64., fmin=0.2, fmax=3.0)
        t = np.arange(0., 2. * np.pi, 360)
        xA, yA = +5., 0.
        xB, yB = -5., 0.
        tS = np.linspace(0., 2.* np.pi, 100)
        xS = 10 * np.cos(tS)
        yS = 10 * np.sin(tS)

        for n in xrange(len(xS)):
            starttime, npts, dt, y = G.ccf(xS[n], yS[n], xA, yA, xB, yB,
                                        delta=0.05,
                                        derivate=True,
                                        geometrical_spreading=True)
            t = starttime + np.arange(npts) * dt
            gca().plot(t, 5. * y + tS[n] * 180./np.pi, "k", alpha=.3)
        gca().set_xlabel('lag time (s)')
        gca().set_ylabel('source azimuth (degrees)')
        showme()

