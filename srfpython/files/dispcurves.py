from srfpython.utils import discrete_time_primitive
from scipy import __version__ as scipyversion
from scipy.interpolate import interp1d
import numpy as np
import copy

assert scipyversion >= "0.14.0"


# ###################################################################################
def C2U(nu, c):
    """convert phase 2 group dispersion curve.
    input :
        nu   = numpy.ndarray, frequency array (positive) of the phase dispersion curve (Hz)
        c    = numpy.ndarray, phase dispersion curve (km/s)
    output :
        nu_  = numpy.ndarray, frequency array of the group dispersion curve (Hz)
        us   = numpy.ndarray, group dispersion curve (km/s)
    """
    Is = nu.argsort()
    nus, cs = nu[Is], c[Is]
    # F = nus / cs
    # us = (nus[1:] - nus[:-1]) / (F[1:] - F[:-1])
    # nu_ = .5 * (nus[1:] + nus[:-1])


    numean = np.sqrt(nus[1:] * nus[:-1])  # see Jeffreys parameters, Tarantola et al 2005
    lu = (np.log(cs[1:]) - np.log(cs[:-1])) / (np.log(nus[1:]) - np.log(nus[:-1]))
    cu = np.interp(numean, xp=nus, fp=cs)
    us = cu / (1. - lu)  # see Lehujet et al., 2017

    return numean, us


# ---------------------------------------
def U2C(nu, u, nuo, co):
    """signaltools.U2C : convert group to phase velocity. Knowledge of the phase dispersion curve is needed at a given frequency nuo.
    nuo must appear into frequency array nu. nuo must be as low as possible.
    input :
        nu   = numpy.ndarray, frequency array (positive) of the group dispersion curve (Hz)
        u    = numpy.ndarray, group dispersion curve (km/s)
        nuo  = float, reference frequency (Hz)
        co   = float, reference phase velocity (km/s)
    output :
        nu_  = numpy.ndarray, frequency array of the phase dispersion curve (Hz)
        c    = numpy.ndarray, phase dispersion curve (km/s)
    """
    Is = nu.argsort()
    nus, us = nu[Is], u[Is]
    if nuo < nu.min() or nuo > nu.max():
        raise ValueError('reference frequency nuo (%f Hz) must be between %f Hz and %f Hz' % (nuo, nu.min(), nu.max()))
    nu_ = nus
    su = discrete_time_primitive(nus, 1. / us, area=False)
    if abs(nu_ - nuo).min() > 0.:
        # need to intepolate su at point nuo
        suo = np.interp(x=nuo, xp=nu_, fp=su)
    else:
        suo = su[abs(nu_ - nuo).argmin()]

    c = co * np.ones(len(nu_))
    Iok = nu_ != nuo
    c[Iok] = nu_[Iok] / (nuo / co + (su[Iok] - suo))
    return nu_, c


# ###################################################################################
class pickable_interp1d(interp1d):
    # I am forced to store the input arguments so that the instance can be pickled and regenerated during depickling
    def __init__(self, x, y, kind='linear', axis=-1, copy=True, bounds_error=True, fill_value=np.nan,
                 assume_sorted=False):
        self.x, self.y = x, y
        self.kind = kind
        self.axis = axis
        self.copy = copy
        self.bounds_error = bounds_error
        self.fill_value = fill_value
        self.assume_sorted = assume_sorted
        self.__setstate__(self.__getstate__())  # trick : __init__ and __setstate__ must keep consistent

    # ---------------------------------------------
    def __getstate__(self):
        return self.x, self.y, self.kind, self.axis, self.copy, self.bounds_error, self.fill_value, self.assume_sorted

    # ---------------------------------------------
    def __setstate__(self, state):
        "By state we denote the tuple of arguments needed to restore the instance dunring depickling"
        x, y, kind, axis, copy, bounds_error, fill_value, assume_sorted = state
        interp1d.__init__(self, x, y, kind=kind, axis=axis, copy=copy, bounds_error=bounds_error,
                          fill_value=fill_value, assume_sorted=assume_sorted)
        # after __init__ because __init__ deletes attributes kind and assumed sorted for unknow reason
        self.kind = kind
        self.assume_sorted = assume_sorted


####################################################################################
class freqinterpolator(object):
    """convention : freqinterpolator(f) = freqinterpolator(-f)"""
    require_positive_values = True
    mode = -1
    flag = "?"
    wave = "?"
    type = "?"
    dvalue = None

    # ---------------------------------------------
    def __init__(self, freq, value, extrapolationmode=1, dvalue=None, **kwargs):
        """
        freq   = numpy.ndarray, frequency values (positive) in Hz
        value  = numpy.ndarray, array corresponding to freq, (postive)
        extrapolationmode : int, 0 1 or 2
            0  = extrapolates by repeating the lowest and highest frequencies = DANGEROUS
            1  = replace extrapoated values by numpy.nan
            2  = raise if the required value is outside the interpolation range
        note  : please do not change attributes manually, use the set_* methods instead
        note2 : freqinterpolator(f) = freqinterpolator(-f)
        """
        assert isinstance(freq, np.ndarray)
        assert isinstance(value, np.ndarray)
        assert len(freq) == len(value)
        assert np.all(freq[1:] > freq[:-1])  # sorted asc
        assert np.all(freq >= 0.)  # positive frequencies only
        if self.require_positive_values:
            assert np.all(value > 0.)

        if dvalue is not None:  # array of uncertaitnies on "value"
            assert isinstance(dvalue, np.ndarray)
            assert len(dvalue) == len(value)
            self.dvalue = dvalue

        self.freq, self.value = freq, value
        self.set_extrapolationmode(extrapolationmode)
        for k, v in kwargs.items(): self.__setattr__(k, v)

    # ---------------------------------------------
    def copy(self):
        return copy.deepcopy(self)

    # ---------------------------------------------
    def set_extrapolationmode(self, extrapolationmode):
        self.extrapolationmode = extrapolationmode  # cannot be changed manually
        if extrapolationmode == 0:
            # this mode is conveniant but dangerous : it extrapolates by repeating values encountered at lowest and highest frequencies,
            # it cannot be used for dispersion curves overtones, because of the cut-off frequencies (see AKI et Richards, 2nd ed, p253)
            freqs = np.concatenate(([-self.freq[0]], self.freq))
            values = np.concatenate(([+self.value[0]], self.value))
            self.i0 = pickable_interp1d(freqs, values, bounds_error=False, fill_value=values[-1], assume_sorted=True)
            self.i1 = self.i2 = None
        elif extrapolationmode == 1:
            # replaces extraplations by np.nan
            self.i1 = pickable_interp1d(self.freq, self.value, bounds_error=False, fill_value=np.nan,
                                        assume_sorted=True)
            self.i0 = self.i2 = None
        elif extrapolationmode == 2:
            # raises in case of extrapolation
            self.i2 = pickable_interp1d(self.freq, self.value, bounds_error=True, fill_value=np.nan, assume_sorted=True)
            self.i0 = self.i1 = None
        else:
            raise ValueError('')

    # ---------------------------------------------
    def call(self, freq):
        if self.extrapolationmode == 1:
            # return np.masked_where(np.isnan(answer), answer)
            return self.i1(abs(freq))
        elif self.extrapolationmode == 2:
            try:
                answer = self.i2(abs(freq))
            except:
                raise Exception('could not evaluate freqinterpolator, check interpolation range')
            return answer
        elif self.extrapolationmode == 0:
            # print "WARNING : extrapolationmode == 0"
            return self.i0(abs(freq))
        else:
            raise Exception('unknown extrapolationmode %d' % self.extrapolationmode)

    # ---------------------------------------------
    def __repr__(self):
        return '%s wave=%s mode=%d type=%s flag=%s extrapmode=%d N=%d' % (
            self.__class__.__name__, self.wave, self.mode, self.type,
            self.flag, self.extrapolationmode, len(self.freq))

    # ---------------------------------------------
    def __call__(self, freq):
        value = self.call(freq)
        if hasattr(freq, "__iter__"): assert len(value) == len(freq)
        return value

    # ---------------------------------------------
    def show(self, ax, marker="-", period=False, showdvalue=True, *args, **kwargs):
        if period:
            hdl = ax.loglog(1. / self.freq, self.value, marker, *args, **kwargs)[0]
        else:
            hdl = ax.loglog(self.freq, self.value, marker, *args, **kwargs)[0]
        if showdvalue and self.dvalue is not None:
            for f, val, dval in zip(self.freq, self.value, self.dvalue):
                # velocity is a Jeffreys parameter
                # val +- dval becomes log(val) +- dval / val
                # then val +- dval becomes exp(log(val) +- dval / val) = val * exp(+- dval / val)
                if "color" not in kwargs.keys():
                    if period:
                        ax.plot([1. / f, 1. / f], [val * np.exp(-dval / val), val * np.exp(+dval / val)], "_-",
                                color=hdl.get_color(), *args, **kwargs)
                    else:
                        ax.plot([f, f], [val * np.exp(-dval / val), val * np.exp(+dval / val)], "_-",
                                color=hdl.get_color(), *args, **kwargs)
                else:
                    if period:
                        ax.plot([1. / f, 1. / f], [val * np.exp(-dval / val), val * np.exp(+dval / val)], "_-", *args,
                                **kwargs)
                    else:
                        ax.plot([f, f], [val * np.exp(-dval / val), val * np.exp(+dval / val)], "_-", *args, **kwargs)

        return hdl


# -------------------------------------------------
class Claw(freqinterpolator):
    """monomodal phase dispersion law
    WARNING : do not use extrapolationmode = 0 for n-th modes because of the cut-off frequencies (see Aki & Richards, 2nd ed, p253)
    """
    require_positive_values = True
    mode = -1  # mode number, 0 = fundamental
    wave = ""  # wave type R, L

    def to_ulaw_old(self):
        gfreq, gvelo = C2U(self.freq, self.value)
        return Ulaw(freq=gfreq, value=gvelo,
                    mode=self.mode, wave=self.wave, extrapolationmode=self.extrapolationmode)

    def to_ulaw(self):
        llaw = self.to_llaw()
        uvalue = self.value / (1. - llaw(self.freq))
        return Ulaw(freq=self.freq, value=uvalue,
                    mode=self.mode, wave=self.wave, extrapolationmode=self.extrapolationmode)


# -------------------------------------------------
class nanClaw(Claw):
    def __init__(self, **kwargs):
        for k, v in kwargs.items(): self.__setattr__(k, v)

    def __call__(self, freq):
        return np.nan * np.ones(len(freq))


# -------------------------------------------------
class Ulaw(Claw):
    "monomodal group dispersion law"
    require_positive_values = True

    def to_claw(self, freq0, c0):
        pfreq, pvelo = U2C(nu=self.freq, u=self.value, nuo=freq0, co=c0)
        return Claw(freq=pfreq, value=pvelo,
                    mode=self.mode, wave=self.wave, extrapolationmode=self.extrapolationmode)


# -------------------------------------------------
class nanUlaw(Ulaw):
    def __init__(self, **kwargs):
        for k, v in kwargs.items(): self.__setattr__(k, v)

    def __call__(self, freq):
        return np.nan * np.ones(len(freq))
