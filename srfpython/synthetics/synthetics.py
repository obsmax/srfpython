from srfpython.depthdisp.dispcurves import Claw, freqspace
from scipy.fftpack import fftfreq, fft, ifft
from signalfuncs import detrend, taperwidth, gaussbandpass, bandpass
import numpy as np

raise Exception('see synthetics2')
# ---------------------------------------
class Green(object):
    """
    simplistic greens function in a 1D elastic dispersive model
    **  account only for the phase dispersion and the gemoetrical spreading of
        a mono modal surface wave
    """
    def __init__(self, claw, fmin=0.2, fmax=5., order=4.0, dt=1. / 12.):
        """
        :param claw: Claw object, phase velocity dispersion laws
        :param fmin: lower frequency for butterworth or central freq for gaussian bandpass, (Hz)
        :param fmax: upper frequency for butterworth or central freq for gaussian bandpass, (Hz)
        :param order: filter order for butterworth or gaussian bandpass
        :param dt: sampling interval (s)
        """
        self.c = claw  #
        self.fmin = fmin  # 0.2
        self.fmax = fmax  # 5.
        self.order = order  # 4.
        self.dt = dt  # 1. / 12.

        assert self.fmax <= 0.5 / self.dt  # nyquist criterion
        if self.fmin < self.fmax:
            f = freqspace(self.fmin, self.fmax, 100, 'flog')
            c = self.c(f)
            self.cmin, self.cmax = np.min(c), np.max(c)
        else:
            self.cmin = self.cmax = self.c(self.fmin)

    # ----------------------------
    def __call__(self, ts, xs, ys, xr, yr):
        """
        compute synthetic waveform
        :param ts: source time, s
        :param xs: source x coordinate, km
        :param ys: source y coordinate, km
        :param xr: receiver x coordinate, km
        :param yr: receiver y coordinate, km
        :return starttime: starttime of the trace in s
        :return npts: number of samples in the trace
        :return dt: sampling interval (s)
        :return y: data array
        """
        distance = ((xs - xr) ** 2. + (ys - yr) ** 2.) ** 0.5
        tmin = np.round(distance / self.cmax - 3. / self.fmin)
        tmax = distance / self.cmin + 3. / self.fmin
        npts = int(np.ceil((tmax - tmin) / self.dt))

        starttime = ts + tmin
        nu = fftfreq(npts, self.dt)
        c = self.c(nu)
        I = ~np.isnan(c)
        Y = np.zeros(len(nu), complex)
        Y[I] = np.exp(-2.j * np.pi * nu[I] * (distance / c[I] - tmin))
        Y /= np.sqrt(distance)  # geometrical spreading
        y = np.real(ifft(Y))

        y = detrend(y)
        y *= taperwidth(y, 1. / self.dt, 3. / self.fmin)
        if self.fmin == self.fmax:
            # gaussian bandpass
            y = gaussbandpass(y, 1. / self.dt, self.fmin, self.order)
        else:
            # bandpass
            y = bandpass(y, 1. / self.dt, self.fmin, self.fmax, self.order, True)

        return starttime, npts, self.dt, y


if __name__ == "__main__":
    from srfpython.standalone.display import *
    from testlaws import c0, c1, c2, c3

    G0 = Green(c0, fmin=0.2, fmax=2., order=4.0, dt=1. / 12.)

    ts, xs, ys = 0., 0., 0.
    for xr in np.arange(1., 20., .25):
        yr = 0.
        starttime, npts, dt, y = G0(ts, xs, ys, xr, yr)
        d = np.sqrt((xs - xr) ** 2. + (ys - yr) ** 2.)
        t = starttime + np.arange(npts) * dt
        gca().plot(t, 4. * y + d, "k", alpha=.3)
    gca().set_xlabel('time (s)')
    gca().set_ylabel('distance (km)')
    showme()