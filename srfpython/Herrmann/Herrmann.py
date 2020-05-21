#!/usr/bin/python2.7
from __future__ import print_function   # python 2/3
from builtins import input   # python 2/3
import os
import signal
import warnings
from subprocess import Popen, PIPE
from numpy import log
import numpy as np
from srfpython.utils import Timer, TimeOutError, Timeout
from srfpython.depthdisp.dispcurves import groupbywtm, igroupbywtm


"""
program to compute surface wave dispersion curves, Maximilien Lehujeur, 01/11/2017
see documentation in function dispersion
use __main__ for demo

WARNING : This module calls fortran programs, make sure they are compiled correctly before calling
"""

_pathfile = os.path.realpath(__file__)  # .../srfpyhon/HerrMann/dispersion.py
_file = _pathfile.split('/')[-1]
_src = _pathfile.rstrip('Herrmann.pyc').rstrip('Herrmann.py') + "src90/"

SRFPRE96_EXE = _pathfile.replace(_file, os.path.join('bin', 'max_srfpre96'))
SRFDIS96_EXE = _pathfile.replace(_file, os.path.join('bin', 'max_srfdis96'))
HERRMANN_TIMEOUT = 5


class CPiSError(Exception):
    pass


class CPiSDomainError(CPiSError):
    pass


class Curve(object):
    def __init__(self, wave, type, mode, freqs, values=None, dvalues=None, checks=True):
        if checks:
            assert wave in ['R', 'L']
            assert type in ['C', 'U']
            assert isinstance(mode, int) and mode >= 0
            assert isinstance(freqs, np.ndarray)
            assert freqs.ndim == 1
            assert np.all(freqs[1:] > freqs[:-1])
        self.wave = wave
        self.type = type
        self.mode = mode
        self.freqs = freqs
        self.values = values
        self.dvalues = dvalues
        self.nfreqs = len(self.freqs)
        self._waves = None  # attribute set only on first demand and kept afterwards
        self._types = None
        self._modes = None

    @property
    def waves(self):
        if self._waves is None:
            self._waves = np.asarray([self.wave for _ in range(self.nfreqs)], '|S1')
        return self._waves

    @property
    def types(self):
        if self._types is None:
            self._types = np.asarray([self.type for _ in range(self.nfreqs)], '|S1')
        return self._types

    @property
    def modes(self):
        if self._modes is None:
            self._modes = np.asarray([self.mode for _ in range(self.nfreqs)], int)
        return self._modes

    def label(self):
        return '{}{}{}'.format(self.wave, self.type, self.mode)

    def plot(self, ax, *args, **kwargs):
        kwargs.setdefault("label", self.label())
        ax.plot(1.0 / self.freqs, self.values, *args, **kwargs)

        ax.set_xlabel('period (s)')
        ax.set_ylabel('velocity (km/s)')
        ax.grid(True, which="major")
        ax.grid(True, which="minor")


def argcrossfind(X, Y):

    """X and Y are unique and sorted"""

    nx, ny = len(X), len(Y)
    ix, iy = 0, 0
    IX, IY = [], []

    while ix < nx and iy < ny:
        if abs((X[ix] - Y[iy]) / X[ix]) < 0.01:
            IX.append(ix)
            IY.append(iy)
            ix += 1
            iy += 1
        elif X[ix] < Y[iy]: ix += 1
        elif X[ix] > Y[iy]: iy += 1
    return np.asarray(IX, int), np.asarray(IY, int)


def readsrfdis96(srfdis96stdout, waves, types, modes, freqs):
    """
    :param srfdis96stdout:
    :param waves:
    :param types:
    :param modes:
    :param freqs:
    :return:
    """
    """converts output from max_srfdis96"""
    periods = 1. / freqs

    # ==== transform the input string
    # remove artifacts whose origin is not known...
    srfdis96stdout = srfdis96stdout.replace(
        '**************', ' 0            ')  # ??? what the fuck

    # strip
    srfdis96stdout = srfdis96stdout.strip().rstrip('\n').split('\n')

    # concatenate rows
    srfdis96stdout = (" ".join(srfdis96stdout)).split()

    # ==== load data
    W   = np.asarray(srfdis96stdout[::6],  int)-1  # wave type 0 = Love, 1 = Rayleigh
    M   = np.asarray(srfdis96stdout[1::6], int)    # (iq-1) mode number, 0 = fundamental
    T1A = np.asarray(srfdis96stdout[2::6], float)  # = t1 if phase only else lower period = t1/(1+h), in s
    T1B = np.asarray(srfdis96stdout[3::6], float)  # = 0. if phase only else upper period = t1/(1-h), in s;
    CC0 = np.asarray(srfdis96stdout[4::6], float)  # = phase velocity at t1 if phase only else at t1a, in km/s;
    CC1 = np.asarray(srfdis96stdout[5::6], float)  # = phase velocity at t1 if phase only else at t1b, in km/s;

    n = len(W)
    I = T1B == 0.  # True means phase only
    L = W == 0     # True means Love
    R = ~L         # assume only R or L

    nI = ~I
    T, C, U = np.zeros(n, float), np.zeros(n, float), np.zeros(n, float) * np.nan
    if I.any():
        # Phase velocity only
        T[I] = T1A[I]
        C[I] = CC0[I]

    if nI.any():
        # Group only or Group and Phase
        T[nI] = 2. * T1A[nI] * T1B[nI] / (T1A[nI] + T1B[nI])
        # see srfpre96 T1A = T1/(1+h), T1B = T1/(1-h) => T1=2 T1A T1B / (T1A + T1B)
        C[nI] = np.sqrt(CC0[nI] * CC1[nI])   # Jeffreys average #A[nI,4:6].mean(axis = 1)
        LnI = (log(CC1[nI]) - log(CC0[nI])) / (log(T1A[nI]) - log(T1B[nI]))
        U[nI] = C[nI] / (1. - LnI)

        # for overtones near the cut-off period
        # one can get same values for CC0 and CC1 => which leads to U=C
        # it seems innacurate
        J = nI & (CC0 == CC1) & (M > 0)
        U[J] = np.nan

    # arange available data
    umodes = np.arange(max(modes) + 1)
    D = {"L": [], "R": []}
    for m in umodes:
        I = (M == m)
        IR = R & I
        IL = L & I
        D["L"].append({"T": T[IL], "C": C[IL], "U": U[IL]})
        D["R"].append({"T": T[IR], "C": C[IR], "U": U[IR]})

    values = np.zeros(len(waves), float) * np.nan
    indexs  = np.arange(len(waves))
    for w, t, m, P, I in groupbywtm(waves, types, modes, periods, indexs, dvalues = None, keepnans = True):
        # available data : period D[w][m]["T"]; value  D[w][m][t]
        # requested data : P
        IP, ITT = argcrossfind(P, D[w][m]["T"])
        values[I[IP]] = D[w][m][t][ITT]

    return values


class HerrmannCallerBasis(object):

    @staticmethod
    def wtmf2srfpre96input(waves, types, modes, freqs):
        assert len(waves) == len(types) == len(modes) == len(freqs)
        fmt = "\nSURF96 {wave:1s} {type:1s} X {mode:d} {period} 1. 1."
        surf96_txt = ""

        nrc = nlc = nru = nlu = 0
        for wave, type, mode, freq in zip(waves, types, modes, freqs):
            surf96_txt += fmt.format(
                wave=wave,
                type=type,
                mode=mode,
                period=1./freq)

            if wave == "R":
                if type == "C":
                    nrc = max([nrc, mode + 1])
                elif type == "U":
                    nru = max([nru, mode + 1])

            elif wave == "L":
                if type == "C":
                    nlc = max([nlc, mode + 1])

                elif type == "U":
                    nlu = max([nlu, mode + 1])

        srfpre96input = """{:d} {:d} {:d} {:d}""".format(nlc, nlu, nrc, nru) \
                        + surf96_txt
        return srfpre96input

    @staticmethod
    def callherrmann(stdin, exe):
        subproc = Popen(
            exe,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE,
            shell=False,
            preexec_fn=os.setsid)

        try:
            with Timeout(HERRMANN_TIMEOUT):
                stdout, stderr = \
                    subproc.communicate(stdin)

        except TimeOutError:
            os.killpg(os.getpgid(subproc.pid), signal.SIGKILL)
            message = "srfpre96 timed out for input :\n" + stdin
            raise CPiSError(message)

        finally:
            try:
                subproc.stdin.close()
                subproc.stdout.close()
                subproc.stderr.close()
            except Exception as e:
                warnings.warn(str(e))
        return stdout, stderr

    @staticmethod
    def depthmodel_arrays_to_string(ztop, vp, vs, rh):
        """prepare input for modified srfpre96 (max_srfpre96)
           z   = depth in km, z[0] must be 0
           vp  = vp in km/s
           vs  = vs in km/s
           rh  = density in g/cm3
        """

        if ztop[0]:
            raise CPiSDomainError('Z0 must be 0')  # assert z[0] == 0.
        n = len(ztop)

        if not (n == len(vp) == len(vs) == len(rh)):
            raise CPiSDomainError('Z VP, VS, RH must be the same length')  # assert n == len(vp) == len(vs) == len(rh)

        ztop = np.asarray(ztop, float)
        if (np.isinf(ztop) | np.isnan(ztop)).any():
            raise CPiSDomainError('got inapropriate values for Z (%s)' % str(ztop))

        Ibad = (ztop[1:] - ztop[:-1] < 0.001)
        if Ibad.any():
            H = ztop[1:] - ztop[:-1]
            raise CPiSDomainError(
                'Z must be growing, layers must be at least '
                '0.001km thick, got %s (%s)' % (str(ztop), str(H[Ibad])))

        vs = np.asarray(vs, float)
        if (np.isinf(vs) | np.isnan(vs)).any():
            raise CPiSDomainError('vs value error %s' % str(vs))

        if not (vs > 0.08).all():
            raise CPiSDomainError('vs domain error %s' % str(vs))

        vp = np.asarray(vp, float)
        if (np.isinf(vp) | np.isnan(vp)).any():
            raise CPiSDomainError('vp value error %s' % str(vp))

        if not (vp > vs * (4. / 3.) ** 0.5).all():
            raise CPiSDomainError('vp over vs domain error %s' % str(vp / vs))

        rh = np.asarray(rh, float)
        if (np.isinf(rh) | np.isnan(rh)).any():
            raise CPiSDomainError('rh value error %s' % str(rh))

        if not np.all((rh > 1.)):
            raise CPiSDomainError('density domain error : %s' % str(rh))

        out = "%d\n" % n + \
              ("%f " * (n - 1)) % tuple(ztop[1:] - ztop[:-1]) + "\n" + \
              ("%f " * n) % tuple(vp) + "\n" + \
              ("%f " * n) % tuple(vs) + "\n" + \
              ("%f " * n) % tuple(rh)
        return out

    def __init__(self, waves, types, modes, freqs, h=0.005, ddc=0.005):
        """
        initiate HerrmannCaller with 4 zippable arrays with same length
        :param waves: 1d array with the wave type letter 'R' for Rayleigh, 'L' for Love
        :param types: 1d array with the wave type letter 'C' for phase vel, 'U' from group vel
        :param modes: 1d array with mode numbers 0=fundamental mode, 1=1st overtone, ...
        :param freqs: 1d array with frequency values, in Hz
        e.g.
        waves = ['R', 'R', 'R', 'L', 'L', 'L']
        types = ['C', 'C', 'C', 'C', 'C', 'C']
        modes = [  0,   0,   0,   0,   0,   0]
        freqs = [ 1.,  2.,  3.,  1.,  2.,  3.]
        :param h: see HerrmannCaller
        :param ddc: see HerrmannCaller
        """

        self.waves = np.asarray(waves, '|S1')
        self.types = np.asarray(types, '|S1')
        self.modes = np.asarray(modes, int)
        self.freqs = np.asarray(freqs, float)

        srfpre96input = self.wtmf2srfpre96input(waves=self.waves, types=self.types, modes=self.modes, freqs=self.freqs)
        self.srfpre96output, stderr = self.callherrmann(stdin=srfpre96input, exe=SRFPRE96_EXE)
        self.srfpre96output = "{:f} {:f}\n".format(h, ddc) + self.srfpre96output

        if len(stderr):
            raise CPiSError('srfpre96 failed with error {}'.format(stderr))

    def disperse(self, ztop, vp, vs, rh):
        """
        :param ztop: top layer depth array, km, first = 0
        :param vp: in km/s, array 1d
        :param vs: in km/s, array 1d
        :param rh: density in g/cm3, array 1d
        :return values: dispersion values in km/s, or nan
        """
        srfdis96input = \
            "{depthmodel_string:s}\n".format(
                depthmodel_string=self.depthmodel_arrays_to_string(ztop, vp, vs, rh)) \
            + self.srfpre96output

        srfdis96output, stderr = self.callherrmann(stdin=srfdis96input, exe=SRFDIS96_EXE)

        if len(stderr):
            raise CPiSError('srfdis96 failed with error : {}'.format(stderr))

        values = readsrfdis96(
            srfdis96stdout=srfdis96output,
            waves=self.waves,
            types=self.types,
            modes=self.modes,
            freqs=self.freqs)
        return values


class HerrmannCaller(HerrmannCallerBasis):
    def __init__(self, curves, h=0.005, ddc=0.005):
        """
        :param curves: a list of Curve objects with the data points at which to compute dispersion
        :param h: factor to convert phase to group see CPS
        :param ddc:
        """

        waves = np.hstack([curve.waves for curve in curves])
        types = np.hstack([curve.types for curve in curves])
        modes = np.hstack([curve.modes for curve in curves])
        freqs = np.hstack([curve.freqs for curve in curves])

        HerrmannCallerBasis.__init__(self, waves=waves, types=types, modes=modes, freqs=freqs, h=h, ddc=ddc)

        self.curves = curves
        self.curve_end_index = np.cumsum([curve.nfreqs for curve in curves])
        self.curve_begin_index = np.concatenate(([0], self.curve_end_index[:-1]))

    def __call__(self, ztop, vp, vs, rh, keepnans=False):
        """
        call Herrmann codes for a depth model
        :param ztop: top layer depth array, km, first = 0
        :param vp: in km/s, array 1d
        :param vs: in km/s, array 1d
        :param rh: density in g/cm3, array 1d
        :param keepnans: whether nans must be kept or not
        :return curves: a list of curves object with the computed values
        """

        values = self.disperse(ztop, vp, vs, rh)

        curves = []
        for b, e in zip(self.curve_begin_index, self.curve_end_index):
            curve = Curve(wave=self.waves[b],
                          type=self.types[b],
                          mode=self.modes[b],
                          freqs=self.freqs[b:e],
                          values=values[b:e])
            curves.append(curve)

        if keepnans:
            return curves

        # remove nans from values
        curves_out = []
        for curve in curves:
            I = np.isnan(curve.values)
            if I.all():
                continue
            curve.freqs, curve.values = curve.freqs[~I], curve.values[~I]
            curves_out.append(curve)
        return curves_out


class HerrmannCallerFromGroupedLists(HerrmannCaller):
    def __init__(self, Waves, Types, Modes, Freqs, h=0.005, ddc=0.005):
        """
        initiate HerrmannCaller with curves data grouped by dispersion curves into arrays
        :param Waves: list of wave letters (one per dispersion curve) ('R' or 'L')
        :param Types: list of wave letters (one per dispersion curve) ('C' or 'U')
        :param Modes: list of mode numbers (one per dispersion curve)
        :param Freqs: list of frequency arrays (one per dispersion curve)
        e.g.
        Waves = ['R', 'L']
        Types = ['C', 'C']
        Modes = [  0,   0]
        Freqs = [ [1.,  2.,  3.],  [1.,  2.,  3.]]
        :param h: see HerrmannCaller
        :param ddc: see HerrmannCaller
        """

        curves = []
        for w, t, m, freqs in zip(Waves, Types, Modes, Freqs):
            curve = Curve(
                wave=w, type=t,
                mode=m, freqs=freqs,
                values=None, dvalues=None)

            curves.append(curve)

        HerrmannCaller.__init__(self, curves=curves, h=h, ddc=ddc)


def check_herrmann_codes():
    """check successful compilation"""

    solution = "please recompile fortran codes using recompile_src90"

    if not os.path.isdir(_src):
        raise Exception('directory %s not found' % _src)

    if not os.path.exists(SRFPRE96_EXE):
        raise Exception('could not find %s\n%s' % (SRFPRE96_EXE, solution))

    if not os.path.exists(SRFDIS96_EXE):
        raise Exception('could not find %s\n%s' % (SRFDIS96_EXE, solution))

    # depth model
    ztop = [0.00, 1.00]  # km, top layer depth
    vp = [2.00, 3.00]  # km/s
    vs = [1.00, 2.00]  # km/s
    rh = [2.67, 2.67]  # g/cm3

    # dipsersion parameters
    curves = [Curve(wave='R', type='C', mode=0, freqs=np.array([1., 2.])),
              Curve(wave='L', type='C', mode=0, freqs=np.array([2., 3.]))]

    hc = HerrmannCaller(curves, h=0.005, ddc=0.005)
    curves = hc(ztop, vp, vs, rh)

    try:
        assert curves[0].wave == "R"
        assert curves[0].type == "C"
        assert curves[0].mode == 0
        assert np.all(curves[0].freqs == np.array([1., 2.]))

        assert curves[1].wave == "L"
        assert curves[1].type == "C"
        assert curves[1].mode == 0
        assert np.all(curves[1].freqs == np.array([2., 3.]))

    except AssertionError:
        raise Exception('could not execute fortran codes\n%s' % solution)


def recompile_src90(yes=False):
    script = """
cd {_src}
./clean.sh && ./compile.sh
""".format(_src=_src)
    if yes:
        os.system(script)
    else:
        print(script)
        if input('run command?').lower() in ["y", "yes"]:
            os.system(script)


def dispersion(*args, **kwargs):
    raise Exception('obsolet')


def dispersion_1(*args, **kwargs):
    raise Exception('obsolet')


def dispersion_2(*args, **kwargs):
    raise Exception('obsolet')


if __name__ == "__main__":

    """ DEMO """

    import time
    import matplotlib.pyplot as plt
    from srfpython.standalone.display import logtick
    check_herrmann_codes()

    # depth model
    ztop = [0.00, 0.25, 0.45, 0.65, 0.85, 1.05, 1.53, 1.80]  # km, top layer depth
    vp   = [1.85, 2.36, 2.63, 3.15, 3.71, 4.54, 5.48, 5.80]  # km/s
    vs   = [0.86, 1.10, 1.24, 1.47, 1.73, 2.13, 3.13, 3.31]  # km/s
    rh   = [2.47, 2.47, 2.47, 2.47, 2.47, 2.58, 2.58, 2.63]  # g/cm3

    # dipsersion parameters
    f = np.logspace(np.log10(0.2), np.log10(3.5), 85)
    curves = [Curve(wave='R', type='U', mode=0, freqs=f),
              Curve(wave='R', type='U', mode=1, freqs=f),
              Curve(wave='R', type='C', mode=0, freqs=f),
              Curve(wave='R', type='C', mode=1, freqs=f),
              Curve(wave='L', type='U', mode=0, freqs=f),
              Curve(wave='L', type='U', mode=1, freqs=f),
              Curve(wave='L', type='C', mode=0, freqs=f),
              Curve(wave='L', type='C', mode=1, freqs=f)]

    hc = HerrmannCaller(curves=curves, h=0.005, ddc=0.005)

    start = time.time()
    curves = hc(ztop=ztop, vp=vp, vs=vs, rh=rh)
    print(time.time() - start)

    # display results
    ax = plt.gca()
    for curve in curves:
        ax.loglog(1. / curve.freqs, curve.values, '+-', label="%s%s%d" % (curve.wave, curve.type, curve.mode))
    ax.set_xlabel('period (s)')
    ax.set_ylabel('velocity (km/s)')
    ax.grid(True, which="major")
    ax.grid(True, which="minor")
    logtick(ax, "xy")
    plt.legend()
    plt.show()
