import numpy as np
from scipy.signal import detrend as scipydetrend
from scipy.fftpack import fft, ifft, fftfreq
from scipy.signal import iirfilter, sosfilt, zpk2sos

def detrend(data):
    return scipydetrend(data, type="linear")


def taperwidth(data, sampling_rate, width):
    tap = np.ones_like(data)
    Nwidth = int(np.round(width * sampling_rate))
    if not Nwidth : return tap

    t = np.arange(len(data)) / sampling_rate
    ttap = 0.5 * (np.sin(np.pi * t[:Nwidth] / float(width) / 1.0 + np.pi / 2.0) + 1.)
    tap[:Nwidth]  *= ttap[::-1]
    tap[-Nwidth:] *= ttap
    return tap


def gaussbandpass(data, sampling_rate, fcenter, alpha):
    npts = len(data)
    dt = 1. / sampling_rate

    nu = fftfreq(npts, dt)
    G = np.exp(-alpha * (np.abs(nu) / fcenter - 1.) ** 2.)
    return np.real(ifft(fft(data) * G))


def bandpass(data, sampling_rate, freqmin, freqmax, corners=4, zerophase=True):
    """modified after obspy
    """
    fe = 0.5 * sampling_rate
    low = freqmin / fe
    high = freqmax / fe
    # raise for some bad scenarios
    if high - 1.0 > -1e-6:
        msg = ("Selected high corner frequency ({}) of bandpass is at or "
               "above Nyquist ({})").format(
            freqmax, fe)
        raise ValueError(msg)
    if low > 1:
        msg = "Selected low corner frequency is above Nyquist."
        raise ValueError(msg)
    z, p, k = iirfilter(corners, [low, high], btype='band',
                        ftype='butter', output='zpk')
    sos = zpk2sos(z, p, k)
    if zerophase:
        firstpass = sosfilt(sos, data)
        return sosfilt(sos, firstpass[::-1])[::-1]
    else:
        return sosfilt(sos, data)
