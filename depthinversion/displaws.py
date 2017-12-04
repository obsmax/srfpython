from utils import pkl, unpkl, munique, discrete_time_primitive, freqspace
from scipy.interpolate import interp1d
import scipy
assert scipy.__version__ >= "0.14.0"

import numpy as np
import copy

####################################################################################
def C2U(nu, c):
    """signaltools.C2U : convert phase 2 group dispersion curve.
    input :
        nu   = numpy.ndarray, frequency array (positive) of the phase dispersion curve (Hz)
        c    = numpy.ndarray, phase dispersion curve (km/s)
    output : 
        nu_  = numpy.ndarray, frequency array of the group dispersion curve (Hz)
        us   = numpy.ndarray, group dispersion curve (km/s)
    """
    Is = nu.argsort()
    nus, cs = nu[Is], c[Is]
    #F = nus / cs
    #us = (nus[1:] - nus[:-1]) / (F[1:] - F[:-1])
    #nu_ = .5 * (nus[1:] + nus[:-1])


    numean = np.sqrt(nus[1:] * nus[:-1]) #see Jeffreys parameters, Tarantola et al 2005
    lu  = (np.log(cs[1:]) - np.log(cs[:-1])) / (np.log(nus[1:]) - np.log(nus[:-1]))
    cu = np.interp(numean, xp = nus, fp = cs)
    us = cu / (1. - lu) #see Lehujet et al., 2017


    return numean, us

####################################################################################
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
    su = discrete_time_primitive(nus, 1./us, area = False) 
    if abs(nu_ - nuo).min() > 0.:  
        #need to intepolate su at point nuo
        suo = np.interp(x = nuo, xp = nu_, fp = su)
    else:
        suo = su[abs(nu_ - nuo).argmin()]

    c = co * np.ones(len(nu_))
    Iok = nu_ != nuo
    c[Iok] = nu_[Iok] / (nuo / co + (su[Iok] - suo))
    return nu_, c

####################################################################################
class pickable_interp1d(interp1d):
    """
    >> old code to be simplified
    """
    def __init__(self, x, y, kind='linear', axis=-1, copy=True, bounds_error=True, fill_value=np.nan, assume_sorted=False):
        self.x, self.y     = x, y
        self.kind          = kind 
        self.axis          = axis
        self.copy          = copy
        self.bounds_error  = bounds_error 
        self.fill_value    = fill_value
        self.assume_sorted = assume_sorted
        self.__setstate__(self.__getstate__()) #trick : __init__ and __setstate__ must keep consistent
    #---------------------------------------------
    def __getstate__(self):
        return (self.x, self.y, self.kind, self.axis, self.copy, self.bounds_error, self.fill_value, self.assume_sorted)
    #---------------------------------------------
    def __setstate__(self, state):
        "By state we denote the tuple of arguments needed to restore the instance during depickling"
        x, y, kind, axis, copy, bounds_error, fill_value, assume_sorted = state
        interp1d.__init__(self, x, y, kind = kind, axis = axis, copy = copy, bounds_error = bounds_error, fill_value = fill_value, assume_sorted = assume_sorted)
        #after __init__ because __init__ deletes attributes kind and assumed sorted for unknow reason
        self.kind          = kind
        self.assume_sorted = assume_sorted

####################################################################################
class freqinterpolator(object):
    """an object to simulate continuous dispersion curve in frequency domain
    convention : freqinterpolator(f) = freqinterpolator(-f)


    >> old code to be simplified
    """
    require_positive_values = True
    mode = -1
    flag = "?"
    wave = "?"
    type = "?"
    dvalue = None
    #---------------------------------------------
    def __init__(self, freq, value, extrapolationmode = 1, dvalue = None, **kwargs):
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
        assert np.all(freq[1:] > freq[:-1])#sorted asc
        assert np.all(freq >= 0.)#positive frequencies only
        if self.require_positive_values:
            assert np.all(value > 0.)


        if dvalue is not None:#array of uncertaitnies on "value"
            assert isinstance(dvalue, np.ndarray)
            assert len(dvalue) == len(value)
            self.dvalue = dvalue

        self.freq, self.value = freq, value
        self.set_extrapolationmode(extrapolationmode)
        for k, v in kwargs.items():self.__setattr__(k, v)

    #---------------------------------------------
    def copy(self):
        return copy.deepcopy(self)

    #---------------------------------------------
    def set_extrapolationmode(self, extrapolationmode):
        self.extrapolationmode = extrapolationmode #cannot be changed manually
        if  extrapolationmode == 0:  
            #this mode is conveniant but dangerous : it extrapolates by repeating values encountered at lowest and highest frequencies, 
            #it cannot be used for dispersion curves overtones, because of the cut-off frequencies (see AKI et Richards, 2nd ed, p253)
            freqs  = np.concatenate(([-self.freq[0]],   self.freq))
            values = np.concatenate(([+self.value[0]], self.value))
            self.i0 = pickable_interp1d(freqs, values, bounds_error = False, fill_value = values[-1], assume_sorted = True)
            self.i1 = self.i2 = None
        elif extrapolationmode == 1: 
            #replaces extraplations by np.nan
            self.i1 = pickable_interp1d(self.freq, self.value, bounds_error = False, fill_value = np.nan, assume_sorted = True)
            self.i0 = self.i2 = None
        elif extrapolationmode == 2: 
            #raises in case of extrapolation
            self.i2 = pickable_interp1d(self.freq, self.value, bounds_error = True,  fill_value = np.nan, assume_sorted = True)        
            self.i0 = self.i1 = None
        else: raise ValueError('')      

    #---------------------------------------------
    def call(self, freq): 
        if self.extrapolationmode == 1: 
            #return np.masked_where(np.isnan(answer), answer)
            return self.i1(abs(freq))
        elif self.extrapolationmode == 2: 
            try: 
                answer = self.i2(abs(freq))
            except:
                raise Exception('could not evaluate freqinterpolator, check interpolation range')
            return answer
        elif self.extrapolationmode == 0 : 
            #print "WARNING : extrapolationmode == 0"
            return self.i0(abs(freq))
        else: raise Exception('unknown extrapolationmode %d' % self.extrapolationmode)

    #---------------------------------------------
    def __repr__(self):
        return '%s wave=%s mode=%d type=%s flag=%s extrapmode=%d N=%d' % (self.__class__.__name__, self.wave, self.mode, self.type, self.flag, self.extrapolationmode, len(self.freq))

    #---------------------------------------------
    def __call__(self, freq):
        value = self.call(freq)
        if hasattr(freq, "__iter__"): assert len(value) == len(freq)
        return value
       
    #---------------------------------------------
    def show(self, ax, marker = "-", period = False, showdvalue = True, *args, **kwargs):
        if period:   hdl = ax.loglog(1./self.freq, self.value, marker, *args, **kwargs)[0]
        else  :      hdl = ax.loglog(   self.freq, self.value, marker, *args, **kwargs)[0]
        if showdvalue and self.dvalue is not None:
            for f, val, dval in zip(self.freq, self.value, self.dvalue):
                if "color" not in kwargs.keys():
                    if period:  ax.plot([1./f, 1./f], [val * np.exp(-dval / val), val * np.exp(+dval / val)], "_-", color = hdl.get_color(), *args, **kwargs)
                    else:       ax.plot([f,       f], [val * np.exp(-dval / val), val * np.exp(+dval / val)], "_-", color = hdl.get_color(), *args, **kwargs)
                else:
                    if period:  ax.plot([1./f, 1./f], [val * np.exp(-dval / val), val * np.exp(+dval / val)], "_-", *args, **kwargs)
                    else:       ax.plot([f,       f], [val * np.exp(-dval / val), val * np.exp(+dval / val)], "_-", *args, **kwargs)  

        return hdl

#-------------------------------------------------
class Claw(freqinterpolator):
    """monomodal phase dispersion law
    WARNING : do not use extrapolationmode = 0 for n-th modes because of the cut-off frequencies (see Aki & Richards, 2nd ed, p253)
    """
    require_positive_values = True
    mode = -1  #mode number, 0 = fundamental
    wave = ""  #wave type R, L

    def to_llaw(self):
        under_freqs = self.freq.copy()
        upper_freqs = self.freq.copy()
        under_freqs[1:]  -= .05 * (under_freqs[1:] - under_freqs[:-1])
        upper_freqs[:-1] += .05 * (upper_freqs[1:] - upper_freqs[:-1])
        lvalue = (np.log10(self(upper_freqs)) - np.log10(self(under_freqs))) / (np.log10(upper_freqs) - np.log10(under_freqs))
        return Llaw(freq = self.freq, value = lvalue, 
            mode = self.mode, wave = self.wave, extrapolationmode = self.extrapolationmode)

    def to_ulaw_old(self):
        gfreq, gvelo = C2U(self.freq, self.value)
        return Ulaw(freq = gfreq, value = gvelo, 
            mode = self.mode, wave = self.wave, extrapolationmode = self.extrapolationmode)

    def to_ulaw(self):
        llaw = self.to_llaw()
        uvalue = self.value / (1. - llaw(self.freq))
        return Ulaw(freq = self.freq, value = uvalue, 
            mode = self.mode, wave = self.wave, extrapolationmode = self.extrapolationmode)    

#-------------------------------------------------
class nanClaw(Claw):
    def __init__(self, **kwargs):
        for k, v in kwargs.items():self.__setattr__(k, v)
    def __call__(self, freq):
        return np.nan * np.ones(len(freq))  

      
#-------------------------------------------------
class Ulaw(Claw):
    "monomodal group dispersion law"
    require_positive_values = True
    def to_claw(self, freq0, c0):
        pfreq, pvelo = U2C(nu = self.freq, u = self.value, nuo = freq0, co = c0)
        return Claw(freq = pfreq, value = pvelo, 
            mode = self.mode, wave = self.wave, extrapolationmode = self.extrapolationmode)

#-------------------------------------------------
class nanUlaw(Ulaw):
    def __init__(self, **kwargs):
        for k, v in kwargs.items():self.__setattr__(k, v)       
    def __call__(self, freq):
        return np.nan * np.ones(len(freq))  

#-------------------------------------------------

class Llaw(freqinterpolator):
    "dlog10(c) / dlog10(f)"
    require_positive_values = False

#----------------------------------------------------
class disppdf(object):
    #2D histogram for dipsersion curves, see also depthmodels.depthpdf, depthmodels.depthpdf2disppdf
    def __init__(self, f, v):
        assert np.all(f[1:] > f[:-1])
        assert np.all(v[1:] > v[:-1])
        self.f = f
        self.v = v
        self.H = np.zeros((len(v), len(f)), float)

    #------------------------------------------------
    def write(self, filename):
        assert filename.split('.')[-1] == "dpdf"
        pkl((self.f, self.v, self.H), filename)

    #------------------------------------------------
    def append(self, law):
        self.appendN(law, Ntimes = 1)

    #------------------------------------------------
    def appendN(self, law, Ntimes=1, **kwargs):
        """append the same model Ntimes times in the histogram"""
        v = law(self.f)
        for i in xrange(len(self.f)):
            if np.isnan(v[i]):continue
            j = np.clip(np.searchsorted(self.v, v[i]), 0, len(self.v) - 1)
            self.H[j, i] += float(Ntimes)
#            if 0 <= j < len(self.v):
#                self.H[j, i] += float(Ntimes)
    #------------------------------------------------
    def appenddat(self, f, v):
        self.appenddatN(f, v, Ntimes = 1)

    #------------------------------------------------
    def appenddatN(self, f, v, Ntimes=1, **kwargs):
        if (f != self.f).any():
            v = np.interp(self.f, xp = f, fp = v, left = np.nan, right = np.nan)
        for i in xrange(len(self.f)):
            if np.isnan(v[i]):continue
            j = np.clip(np.searchsorted(self.v, v[i]), 0, len(self.v) - 1)
            self.H[j, i] += float(Ntimes)
    #------------------------------------------------
    def show(self, ax, **kwargs):
        I = np.zeros(len(self.f))
        for i in xrange(len(self.f)):
            I[i] = discrete_time_primitive(self.v, self.H[:, i], area = True)
        I = discrete_time_primitive(self.f, I, area = True)
        H = np.ma.masked_where(self.H == 0., self.H)
        ax.pcolormesh1(1. / self.f, self.v, H / I, **kwargs)   
        ax.set_xscale('log')
        ax.set_yscale('log')
    #------------------------------------------------
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
#------------------------------------------------
class disppdf_from_zpdffile(disppdf):
    def __init__(self, filename):
        self.f, self.v, self.H = unpkl(filename)

#---------------------------------------------
class surf96reader_from_s96_txt(object):
    """>> to be simplified"""
    def __init__(self, s96txt):

        self.data = {"WAVE": [], "TYPE": [], "FLAG": [], "MODE": [], "PERIOD": [], "VALUE": [], "DVALUE": []}
        self.filename = ""
        for l in s96txt.split('\n'):
            if not len(l.strip()):continue
            if not l.startswith('SURF96'):raise Exception('line did not start with SURF96')
            l = l.split()
            if not len(l) == 8: raise Exception('incorrect number of fields in line')
            _, wave, typ, flag, mode, per, value, dvalue = l
            mode = int(mode)
            per, value, dvalue = map(float, [per, value, dvalue])

            self.data['WAVE'].append(wave)
            self.data['TYPE'].append(typ)
            self.data['FLAG'].append(flag)
            self.data['MODE'].append(mode)
            self.data['PERIOD'].append(per)
            self.data['VALUE'].append(value)
            self.data['DVALUE'].append(dvalue)

        for k in self.data.keys():
            self.data[k] = np.asarray(self.data[k])

        
    #---------------------------------------------
    def get_all(self):
        for w, t, f, m in zip(*munique(self.data['WAVE'], self.data['TYPE'], self.data['FLAG'], self.data['MODE'])):
            yield self.get(mode=m, wave=w, type=t, flag=f)


    #---------------------------------------------
    def get(self, wave = "R", type = "C", mode = 0,  flag = None):
        """generate a dispersion law according to the (read) file content"""
        I = (self.data['MODE'] == mode) & (self.data['WAVE'] == wave) & (self.data['TYPE'] == type)
        if flag is not None: I = I & (self.data['FLAG'] == flag)

        if not np.any(I): 
            raise Exception('mode %d, wave %s, type %s and flag %s not found in %s' % (mode, wave, type, flag, self.filename))
        freq = 1. / self.data['PERIOD'][I]
        valu = self.data['VALUE'][I]
        dval = self.data['DVALUE'][I]
        J = np.argsort(freq)       
        if   type == "C": W = Claw
        elif type == "U": W = Ulaw
        else: raise ValueError('only for types C and U')
        return W(freq[J], valu[J], extrapolationmode = int(mode != 0), dvalue = dval[J], mode = mode, wave = wave, type = type, flag = flag)

    #---------------------------------------------
    def get_lower(self, wave = "R", type = "C", mode = 0,  flag = None):
        law = self.get(wave = wave, type = type, mode = mode,  flag = flag)
        if   type == "C": W = Claw
        elif type == "U": W = Ulaw
        else: raise ValueError('only for types C and U')
        return W(law.freq, law.value * np.exp(-law.dvalue / law.value), extrapolationmode = int(mode != 0), dvalue = None, mode = mode, wave = wave, type = type, flag = flag)

    #---------------------------------------------
    def get_upper(self, wave = "R", type = "C", mode = 0,  flag = None):
        law = self.get(wave = wave, type = type, mode = mode,  flag = flag)
        if   type == "C": W = Claw
        elif type == "U": W = Ulaw
        else: raise ValueError('only for types C and U')
        return W(law.freq, law.value * np.exp(+law.dvalue / law.value), extrapolationmode = int(mode != 0), dvalue = None, mode = mode, wave = wave, type = type, flag = flag)

    #---------------------------------------------
    def wtm(self):
        for w, t, m in zip(*munique(self.data['WAVE'], self.data['TYPE'], self.data['MODE'])):
            yield w, t, m
    #---------------------------------------------
    def wtmfvd(self):
        return self.data['WAVE'], \
               self.data['TYPE'], \
               self.data['MODE'], \
               1. / self.data['PERIOD'], \
               self.data['VALUE'], \
               self.data['DVALUE']
    #---------------------------------------------
    def wtmfv(self):
        return self.data['WAVE'], \
               self.data['TYPE'], \
               self.data['MODE'], \
               1. / self.data['PERIOD'], \
               self.data['VALUE']
    #---------------------------------------------
    def __str__(self):
        fmt = "SURF96 {wave} {type} {flag} {mode} {period} {value} {dvalue}"
        return "\n".join([fmt.format(\
                wave = w,
                type = t,
                mode = m,
                period = p,
                value = v,
                dvalue = d,
                flag = "X") 
            for w, t, m, p, v, d in \
                zip(self.data['WAVE'],
                    self.data['TYPE'],
                    self.data['MODE'],
                    self.data['PERIOD'],
                    self.data['VALUE'],
                    self.data['DVALUE'])])
    # ---------------------------------------------------
    def write96(self, filename):
        with open(filename, "w") as fid:
            fid.write(self.__str__())

    # ---------------------------------------------
    def copy(self):
        return copy.deepcopy(self)
        
#-------------------------------------------------
class surf96reader(surf96reader_from_s96_txt):
    def __init__(self, filename):
        with open(filename, 'r') as fid:
            surf96reader_from_s96_txt.__init__(self, "".join(fid.readlines()))
        self.filename = filename



#-------------------------------------------------
if __name__ == "__main__":
    from labex.graphictools.gutils import *
    txt = """SURF96 R U X   0 20.000000 2.980498 0.298050
SURF96 R U X   0 15.695199 2.907702 0.290770
SURF96 R U X   0 12.316964 2.808471 0.280847
SURF96 R U X   0 9.665860 2.667765 0.266776
SURF96 R U X   0 7.585380 2.458498 0.245850
SURF96 R U X   0 5.952703 2.139880 0.213988
SURF96 R U X   0 4.671443 1.741809 0.174181
SURF96 R U X   0 3.665961 1.476472 0.147647
SURF96 R U X   0 2.876900 0.488186 0.048819
SURF96 R U X   0 2.257676 0.527731 0.052773
SURF96 R U X   0 1.771734 0.744356 0.074436
SURF96 R U X   0 1.390386 0.872462 0.087246
SURF96 R U X   0 1.091119 0.924153 0.092415
SURF96 R U X   0 0.856266 0.885369 0.088537
SURF96 R U X   0 0.671964 0.607388 0.060739
SURF96 R U X   0 0.527330 0.328387 0.032839
SURF96 R U X   0 0.413828 0.337383 0.033738
SURF96 R U X   0 0.324755 0.398124 0.039812
SURF96 R U X   0 0.254855 0.438161 0.043816
SURF96 R U X   0 0.200000 0.458669 0.045867
SURF96 R U X   1 4.220795 1.707556 0.170756
SURF96 R U X   1 3.563021 1.245423 0.124542
SURF96 R U X   1 3.007756 0.934693 0.093469
SURF96 R U X   1 2.539024 1.422822 0.142282
SURF96 R U X   1 2.143340 1.433326 0.143333
SURF96 R U X   1 1.809320 1.413521 0.141352
SURF96 R U X   1 1.527353 1.367689 0.136769
SURF96 R U X   1 1.289329 1.271752 0.127175
SURF96 R U X   1 1.088398 1.075943 0.107594
SURF96 R U X   1 0.918781 0.843961 0.084396
SURF96 R U X   1 0.775597 0.675538 0.067554
SURF96 R U X   1 0.654727 0.627287 0.062729
SURF96 R U X   1 0.552694 0.754907 0.075491
SURF96 R U X   1 0.466562 0.746757 0.074676
SURF96 R U X   1 0.393852 0.698782 0.069878
SURF96 R U X   1 0.332474 0.661084 0.066108
SURF96 R U X   1 0.280661 0.631484 0.063148
SURF96 R U X   1 0.236922 0.585542 0.058554
SURF96 R U X   1 0.200000 0.472077 0.047208
SURF96 R C X   0 20.000000 3.100898 0.310090
SURF96 R C X   0 15.695199 3.065806 0.306581
SURF96 R C X   0 12.316964 3.018683 0.301868
SURF96 R C X   0 9.665860 2.953747 0.295375
SURF96 R C X   0 7.585380 2.860367 0.286037
SURF96 R C X   0 5.952703 2.717703 0.271770
SURF96 R C X   0 4.671443 2.496378 0.249638
SURF96 R C X   0 3.665961 2.221357 0.222136
SURF96 R C X   0 2.876900 1.743041 0.174304
SURF96 R C X   0 2.257676 1.078402 0.107840
SURF96 R C X   0 1.771734 0.940179 0.094018
SURF96 R C X   0 1.390386 0.910435 0.091043
SURF96 R C X   0 1.091119 0.909286 0.090929
SURF96 R C X   0 0.856266 0.910326 0.091033
SURF96 R C X   0 0.671964 0.873889 0.087389
SURF96 R C X   0 0.527330 0.704534 0.070453
SURF96 R C X   0 0.413828 0.560508 0.056051
SURF96 R C X   0 0.324755 0.504186 0.050419
SURF96 R C X   0 0.254855 0.483503 0.048350
SURF96 R C X   0 0.200000 0.475880 0.047588
SURF96 R C X   1 4.220795 3.351892 0.335189
SURF96 R C X   1 3.563021 2.806008 0.280601
SURF96 R C X   1 3.007756 2.133239 0.213324
SURF96 R C X   1 2.539024 1.936461 0.193646
SURF96 R C X   1 2.143340 1.835743 0.183574
SURF96 R C X   1 1.809320 1.756824 0.175682
SURF96 R C X   1 1.527353 1.688001 0.168800
SURF96 R C X   1 1.289329 1.618637 0.161864
SURF96 R C X   1 1.088398 1.529330 0.152933
SURF96 R C X   1 0.918781 1.395882 0.139588
SURF96 R C X   1 0.775597 1.230463 0.123046
SURF96 R C X   1 0.654727 1.071948 0.107195
SURF96 R C X   1 0.552694 0.987828 0.098783
SURF96 R C X   1 0.466562 0.943759 0.094376
SURF96 R C X   1 0.393852 0.900495 0.090049
SURF96 R C X   1 0.332474 0.856767 0.085677
SURF96 R C X   1 0.280661 0.815215 0.081521
SURF96 R C X   1 0.236922 0.774839 0.077484
SURF96 R C X   1 0.200000 0.723221 0.072322"""
    s = surf96reader_from_s96_txt(txt)
    for _ in s.get_all():
        print _
        _.show(gca(), period=True)
    showme()


#-------------------------------------------------
# def premklaws(waves, types, modes, freqs, values, dvalues = None, keepnans = False):
#     """group dispersion data by wave, type, mode
#     """
#     freqs = np.array(freqs, float)
#     types = np.array(list(types), "|S1")
#     waves = np.array(list(waves), "|S1")
#     modes = np.array(modes,  int)
#     assert len(waves) == len(types) == len(modes) == len(freqs)
#     if dvalues is not None:
#         dvalues = np.array(dvalues, float)
#         assert len(waves) == len(dvalues)
#
#     w, t, m = munique(waves, types, modes)
#     I = np.lexsort((m, t, w))
#
#     for w, t, m in zip(w[I], t[I], m[I]):
#         J = (waves == w) & (types == t) & (modes == m)
#
#         if keepnans:   K = np.ones(len(values[J]), bool)
#         else:          K = ~np.isnan(values[J])
#
#         L = np.argsort(freqs[J][K])
#
#         if dvalues is None:
#             yield w, t, m, freqs[J][K][L], values[J][K][L]
#         else:
#             yield  w, t, m, freqs[J][K][L], values[J][K][L], dvalues[J][K][L]
# #-------------------------------------------------
# def mklaws(waves, types, modes, freqs, values, dvalues = None):
#     """values are from disperse"""
#     laws = []
#     for tup in premklaws(waves, types, modes, freqs, values, dvalues = dvalues, keepnans = False):
#         if dvalues is None:
#             w, t, m, F, V = tup
#             D = None
#         else :
#             w, t, m, F, V, D = tup
#
#         whichlaw = [[nanClaw, nanUlaw], [Claw, Ulaw]][len(F) > 1][t == "U"]
#         law = whichlaw(freq = F, value = V, dvalue = D, extrapolationmode = int(m > 0), mode = m, type = t, wave = w)
#         laws.append(law)
#     return laws

##-------------------------------------------------
#class sdisprascreader(Allfile):
#    """initiate dispersion law from Hermann's file SDISPR.ASC : contains both phase speed and energy!!"""
#    def __init__(self, sdisprasc):
#        self.filename = os.path.realpath(sdisprasc)
#        wave = "R"#see SLEGN.ASC for love waves
#        if not os.path.exists(sdisprasc): raise FileNotFound('%s not found' % sdisprasc)
#        with open(sdisprasc, 'r') as fid:
#            L = fid.readlines()
#            if not len(L): raise EmptyFile("%s is empty" % sdisprasc,)
#            L = L[1:]#skip header line
#            mode, nfreq, period, frequency, phasespeed = np.zeros(len(L), int), np.zeros(len(L), int), np.zeros(len(L), float), np.zeros(len(L), float), np.zeros(len(L), float)
#            for i, l in enumerate(L):
#                m, nf, p, f, c = l.split("\n")[0].split()
#                mode[i], nfreq[i], period[i], frequency[i], phasespeed[i] = int(m), int(nf), abs(float(p)), abs(float(f)), abs(float(c))

#        Allfile.__init__(self, '')
#        self.data.MODE  = mode;  self.fmt.MODE  = "%d"
#        self.data.FREQ  = frequency;  self.fmt.FREQ  = "%f"
#        self.data.PHASE = phasespeed; self.fmt.PHASE = "%f"

#    #---------------------------------------------
#    def get(self, mode, type = "C"):
#        I = self.data.MODE == mode
#        if not sum(I): 
#            os.system('cat %s' % self.filename)
#            raise Exception('mode %d and type %s not found if %s' % (mode, type, self.filename))
#        freq = self.data.FREQ[I]
#        valu = self.data.PHASE[I]
#        J = np.argsort(freq)   
#        claw = Claw(freq[J], valu[J], extrapolationmode = int(mode != 0), mode = mode, wave = "R", type = "C")
#        if type == "C":
#            return claw
#        elif type == "U":
#            try:    ulaw = claw.to_ulaw()
#            except: raise CPiSError('')
#            return ulaw
#        else: raise
##-------------------------------------------------

#if __name__ == "__main__":
#    from labex.graphictools.gutils import *
#    from labex.tools.generic_nobspy import unpkl
#    #s = surf96reader('/home/max/progdat/CPiS/EarthModel/Soultz.surf96')
#    #c = s.get()
#    #c.show(gca())
#    #showme()
#    c = unpkl('/home/max/progdat/CPiS/EarthModel/clawSZRC0.pkl')

#    c.show(gca())
#    showme()    


