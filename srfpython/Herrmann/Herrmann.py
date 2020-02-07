#!/usr/bin/python2.7
from __future__ import print_function
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

SRFPRE96_EXE = _pathfile.replace(_file, 'bin/max_srfpre96')
SRFDIS96_EXE = _pathfile.replace(_file, 'bin/max_srfdis96')
#SHELL_COMMAND = "{}|{}".format(SRFPRE96_EXE, SRFDIS96_EXE)
HERRMANN_TIMEOUT = 5

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
    curves = [('R', 'C', 0, np.array([1., 2.])),
              ('L', 'C', 0, np.array([2., 3.]))]

    g = dispersion_2(ztop, vp, vs, rh, curves)

    try:
        w, t, m, F, V = g.next()
        assert (w, t, m) == ("L", "C", 0)
        assert np.all(F == np.array([2., 3.]))
        w, t, m, F, V = g.next()
        assert (w, t, m) == ("R", "C", 0)
        assert np.all(F == np.array([1., 2.]))

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
        if raw_input('run command?').lower() in ["y", "yes"]:
            os.system(script)


# ##################################### MODIFIED HERRMANN'S CODES

def prep_srfpre96_1(h=0.005, dcl=0.005, dcr=0.005):
    """prepare input for modified srfpre96 (max_srfpre96)
    dcl, dcr are root search increment for loev and rayleigh waves respectively
    h, dcl,dcr : see srfpre96
    """
    return "%f %f %f" % (h, dcl, dcr)


def prep_srfpre96_2(z, vp, vs, rh):
    """prepare input for modified srfpre96 (max_srfpre96)
       z   = depth in km, z[0] must be 0
       vp  = vp in km/s
       vs  = vs in km/s
       rh  = density in g/cm3
    """

    if z[0]:
        raise CPiSDomainError('Z0 must be 0')  # assert z[0] == 0.
    n = len(z)

    if not (n == len(vp) == len(vs) == len(rh)):
        raise CPiSDomainError('Z VP, VS, RH must be the same length')  # assert n == len(vp) == len(vs) == len(rh)

    z = np.asarray(z, float)
    if (np.isinf(z) | np.isnan(z)).any():
        raise CPiSDomainError('got inapropriate values for Z (%s)' % str(z))

    Ibad = (z[1:] - z[:-1] < 0.001)
    if Ibad.any():
        H = z[1:] - z[:-1]
        raise CPiSDomainError(
            'Z must be growing, layers must be at least '
            '0.001km thick, got %s (%s)' % (str(z), str(H[Ibad])))

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
          ("%f " * (n - 1)) % tuple(z[1:] - z[:-1]) + "\n" + \
          ("%f " * n) % tuple(vp) + "\n" + \
          ("%f " * n) % tuple(vs) + "\n" + \
          ("%f " * n) % tuple(rh)
    return out


def prep_srfpre96_3(waves,types,modes,freqs):
    """prepare input for modified srfpre96 (max_srfpre96)"""
    nlc = nlu = nrc = nru = 0
    fmt="\nSURF96 %1s %1s X %d %f 1. 1."
    out = ""

    for w,t,m,f in zip(waves,types,modes,freqs):
        out += fmt % (w,t,m,1./f)

        if w == "R":
            if t == "C":
                nrc = max([nrc, m + 1])
            elif t == "U":
                nru = max([nru, m + 1])

        elif w == "L":
            if t == "C":
                nlc = max([nlc, m + 1])

            elif t == "U":
                nlu = max([nlu, m + 1])

    out = "%d %d %d %d" %(nlc,nlu,nrc,nru) + out
    return out


class CPiSError(Exception):
    pass


class CPiSDomainError(CPiSError):
    pass


def readsrfdis96_old_stable(stdout, waves, types, modes, freqs):
    """converts output from max_srfdis96, to be optimized
    (takes about 0.5 * the times used by max_srfpre96+max_srfdis96)"""

    # try:
    #     A = np.genfromtxt(StringIO(unicode(stdout)), dtype = float)
    # except:
    #     raise CPiSError('could not understand stdout \n%s' % stdout)

    A = np.asarray((" ".join(stdout.strip().rstrip('\n').split('\n'))).split(), float)
    A = A.reshape((len(A) / 6, 6))

    # A[:, 0] (itst) wave type 1 = Love, 2 = Rayleigh
    # A[:, 1] (iq-1) mode number, 0 = fundamental
    # A[:, 2] (t1a)  t1 if phase only else lower period = t1/(1+h), in s
    # A[:, 3] (t1b)  0. if phase only else upper period = t1/(1-h), in s;
    # A[:, 4] (cc0)  phase velocity at t1 if phase only else at t1a, in km/s;
    # A[:, 5] (cc1)  phase velocity at t1 if phase only else at t1b, in km/s;

    W = A[:, 0]
    M = A[:, 1]
    I = A[:, 3] == 0.  # True means phase only
    nI = ~I
    n = A.shape[0]
    T, C, U = np.zeros(n, float), np.zeros(n, float), np.zeros(n, float) * np.nan

    if I.any():
        T[I]  = A[I,2]
        C[I]  = A[I,4]

    if nI.any():
        # np.sqrt(A[nI,2] * A[nI,3]) #Jeffreys average #A[nI,2:4].mean(axis = 1)
        T[nI] = A[nI, 2] * A[nI, 3] / (A[nI, 2:4].mean(axis=1))

        C[nI] = np.sqrt(A[nI,4] * A[nI,5])  # Jeffreys average  # A[nI,4:6].mean(axis = 1)

        LnI = (log(A[nI,5]) - log(A[nI,4])) / (log(A[nI,2]) - log(A[nI,3]))
        U[nI] = C[nI] / (1. - LnI)

    L, R = (W == 1), (W == 2)
    umodes = np.arange(max(modes) + 1)
    RMs = [R & (M == m) for m in umodes]
    LMs = [L & (M == m) for m in umodes]
    RMs = [rms if rms.any() else None for rms in RMs]
    LMs = [lms if lms.any() else None for lms in LMs]

    values = np.zeros(len(waves), float) * np.nan
    for n, (w, t, m, f) in enumerate(zip(waves, types, modes, freqs)):

        if w == "R":
            S = RMs[m]
        else:
            S = LMs[m]  # elif w == "L"
        if S is None:
            continue

        p  = 1. / f
        TS = T[S]
        iS = np.abs(TS - p).argmin()
        per = TS[iS]
        if abs(per - p) / p > 0.01: continue

        if t == 'C':
            val = C[S][iS]

        else:
            val = U[S][iS] #elif t == "U"
        # else:
        #      raise ValueError('')
        # if val <= 0: #NON, je met une penalite dans la fonction cout
        #     raise CPiSError('encountered negative dispersion velocity')

        values[n] = val
    return values


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


def readsrfdis96(stdout, waves, types, modes, freqs):
    """converts output from max_srfdis96"""
    periods = 1./freqs

    stdout = stdout.replace('**************', ' 0            ') #??? what the fuck
    stdout = stdout.strip().rstrip('\n').split('\n')
    stdout = [_[:2] + " " + _[2:] for _ in stdout] # add one more space for mode numbers higher than 10

    A = (" ".join(stdout)).split()
    W   = np.asarray(A[::6],  int)-1  # wave type 0 = Love, 1 = Rayleigh
    M   = np.asarray(A[1::6], int)    # (iq-1) mode number, 0 = fundamental
    T1A = np.asarray(A[2::6], float)  # t1 if phase only else lower period = t1/(1+h), in s
    T1B = np.asarray(A[3::6], float)  # 0. if phase only else upper period = t1/(1-h), in s;
    CC0 = np.asarray(A[4::6], float)  # phase velocity at t1 if phase only else at t1a, in km/s;
    CC1 = np.asarray(A[5::6], float)  # phase velocity at t1 if phase only else at t1b, in km/s;

    n = len(W)
    I = T1B == 0. #True means phase only
    L = W == 0    #True means Love
    R = ~L        #assume only R or L

    nI = ~I
    T, C, U = np.zeros(n, float), np.zeros(n, float), np.zeros(n, float) * np.nan
    if I.any():
        T[I] = T1A[I]
        C[I] = CC0[I]

    if nI.any():
        T[nI] = 2. * T1A[nI] * T1B[nI] / (T1A[nI] + T1B[nI])
        # see srfpre96 T1A = T1/(1+h), T1B = T1/(1-h) => T1=2 T1A T1B / (T1A + T1B)
        C[nI] = np.sqrt(CC0[nI] * CC1[nI])   # Jeffreys average #A[nI,4:6].mean(axis = 1)
        LnI = (log(CC1[nI]) - log(CC0[nI])) / (log(T1A[nI]) - log(T1B[nI]))
        U[nI] = C[nI] / (1. - LnI)

    # arange available data
    umodes = np.arange(max(modes) + 1)
    D = {"L": [], "R": []}
    for m in umodes:
        I = (M == m)
        IR = R & I
        IL = L & I
        D["L"].append({"T" : T[IL], "C" : C[IL], "U" : U[IL]})
        D["R"].append({"T" : T[IR], "C" : C[IR], "U" : U[IR]})

    values = np.zeros(len(waves), float) * np.nan
    indexs  = np.arange(len(waves))
    for w, t, m, P, I in groupbywtm(waves, types, modes, periods, indexs, dvalues = None, keepnans = True):
        # available data : period D[w][m]["T"]; value  D[w][m][t]
        # requested data : P
        IP, ITT = argcrossfind(P, D[w][m]["T"])
        values[I[IP]] = D[w][m][t][ITT]

    return values


def dispersion(ztop, vp, vs, rh,
    waves, types, modes, freqs,
    h=0.005, dcl=0.005, dcr=0.005):

    """compute surface wave dispersion curves from a 1-D depth model
    based on a modified version of the codes from Herrmann and Ammon 2002

    *) a dispersion curve is given by 4 attributes
        wave : string, "R" for Rayleigh, "L" for Love
        type : string, "C" for Phase velocity, "U" for Group-velocity
        mode : integer, a mode number, 0 means fundamental
        freq : array, frequencies in Hz

    *) a 1-D depth model is given by 4 attributes
        ztop : list or array, sorted, positive, top layer depth in km, ztop[0] must be 0 !!!
        vp   : list or array, P wave velocity in km/s
        vs   : list or array, S wave velocity in km/s
        rh   : list or array, density in g.cm^-3

    input:
        -> depth model
        ztop, vp, vs, rh = depth model, 4 iterables with same length

        -> dispersion points
        waves, types, modes, freqs = dispersion curves, 4 iterables with same length
        example :
            waves = ['L', 'L', 'L', 'R', 'R', 'R', 'R', 'R']
            types = ['U', 'U', 'U', 'C', 'C', 'C', 'C', 'C']
            modes = [ 0,   0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ]
            freqs = [ 1.,  2.,  3.,  1.,  2.,  3.,  4.,  5.]

        -> Herrmann's parameters, see CPS documentation
        h = float, period increment for phase to group conversion (0.005 is reasonable)
        dcl, dcr = 2 floats, phase velocity increment for love and rayleigh root searches respectively, see Herrmann doc

    output:
        values : list or array, dispersion values for each wave, type, mode, possibly nan (above cut-off period)
                 note : use groupbywtm to group outputs by wave, type, and mode
    see also :
        dispersion_1
        dispersion_2
        groupbywtm
        igroupbywtm
    """

    instr2 = prep_srfpre96_2(ztop, vp, vs, rh)
    instr1 = prep_srfpre96_1(h=h, dcl=dcl, dcr=dcr)
    instr3 = prep_srfpre96_3(waves, types, modes, freqs)
    pstdin = "\n".join([instr2, instr1, instr3])

    SRFPRE96_SUBPROC = Popen(SRFPRE96_EXE,
                             stdin=PIPE,
                             stdout=PIPE,
                             stderr=PIPE,
                             shell=False,
                             preexec_fn=os.setsid)
    SRFDIS96_SUBPROC = Popen(SRFDIS96_EXE,
                             stdin=PIPE,  # SRFPRE96_SUBPROC.stdout,
                             stdout=PIPE,
                             stderr=PIPE,
                             shell=False,
                             preexec_fn=os.setsid)
    try:
        with Timeout(HERRMANN_TIMEOUT):
            pstdout, pstderr = SRFPRE96_SUBPROC.communicate(pstdin)
            qstdout, qstderr = SRFDIS96_SUBPROC.communicate(pstdout)

    except TimeOutError as e:
        os.killpg(os.getpgid(SRFPRE96_SUBPROC.pid), signal.SIGKILL)
        os.killpg(os.getpgid(SRFDIS96_SUBPROC.pid), signal.SIGKILL)
        message = "Herrmann timed out for model:\n\tztop={}\n\tvp={}\n\tvs={}\n\trh={}\n".format(ztop, vp, vs, rh)
        raise CPiSError(message)

    finally:
        try:
            SRFPRE96_SUBPROC.stdin.close()
            SRFDIS96_SUBPROC.stdin.close()
            SRFPRE96_SUBPROC.stdout.close()
            SRFDIS96_SUBPROC.stdout.close()
            SRFPRE96_SUBPROC.stderr.close()
            SRFDIS96_SUBPROC.stderr.close()
        except Exception as e:
            warnings.warn(str(e))

    values = readsrfdis96(qstdout, waves, types, modes, freqs)
    return values


def dispersion_1(ztop, vp, vs, rh,
    Waves, Types, Modes, Freqs,
    h = 0.005, dcl = 0.005, dcr = 0.005, keepnans = False):

    """same as dispersion with slightly more convenient inputs and outputs
    (no need to repeat wave, type and mode)

    input:
        -> depth model
        ztop, vp, vs, rh = depth model, 4 iterables with same length

        -> dispersion points
        Waves, Types, Modes, Freqs = dispersion curves, 4 iterables with same length
        example :
            Waves = ['L', 'L', 'R']
            Types = ['C', 'C', 'U']
            Modes = [ 0,   1,   0 ]
            Freqs = [fLC0, fLC1, fRU0] where f??? are frequency numpy arrays or lists

        -> Herrmann's parameters, see dispersion

    output :
        a generator of tuples
        each tuple correspond to one dispersion curve
        (wave(str,"R"/"L"), type(str,"C"/"U"), mode(int), frequency(array,Hz), velocity(array,km/s))
        example
            ("L", "C", 0, fLC0, VLC0)
            ("L", "C", 1, fLC1, VLC1)
            ("R", "U", 0, fRU0, VRU0)

    see also :
        dispersion
        dispersion_2
        groupbywtm
        igroupbywtm
    """
    waves, types, modes, freqs = igroupbywtm(Waves, Types, Modes, Freqs)
    values = dispersion(ztop, vp, vs, rh, \
                waves, types, modes, freqs,
                h = h, dcl = dcl, dcr = dcr)

    for w, t, m, F, V in groupbywtm(waves, types, modes, freqs, values, keepnans = keepnans):
        yield w, t, m, F, V


def dispersion_2(ztop, vp, vs, rh, Curves,
    h=0.005, dcl=0.005, dcr=0.005, keepnans=False):

    """same as dispersion with slightly more convenient inputs and outputs
    (inputs are grouped by dispersion curves)

    input:
        -> depth model
        ztop, vp, vs, rh = depth model, 4 iterables with same length

        -> dispersion curves
        Curves = list of tuples like (wave(str,"R"/"L"), type(str,"C"/"U"), mode(int), frequency(array,Hz))
        example :
            Curves = [('L', 'C', 0, fLC0),
                      ('L', 'C', 1, fLC1),
                      ('R', 'U', 0, fRU0)]

        -> Herrmann's parameters, see dispersion

    output :
        see dispersion_1

    see also :
        dispersion
        dispersion_1
        groupbywtm
        igroupbywtm
    """
    Waves, Types, Modes, Freqs = zip(*Curves)
    for w, t, m, F, V in dispersion_1(\
            ztop, vp, vs, rh,
            Waves, Types, Modes, Freqs,
            h=h, dcl=dcl, dcr=dcr, keepnans=keepnans):
        yield w, t, m, F, V


if __name__ == "__main__":

    """ DEMO """

    import matplotlib.pyplot as plt
    from srfpython.standalone.display import logtick
    check_herrmann_codes()

    # depth model
    ztop = [0.00, 0.25, 0.45, 0.65, 0.85, 1.05, 1.53, 1.80]  # km, top layer depth
    vp   = [1.85, 2.36, 2.63, 3.15, 3.71, 4.54, 5.48, 5.80]  # km/s
    vs   = [0.86, 1.10, 1.24, 1.47, 1.73, 2.13, 3.13, 3.31]  # km/s
    rh   = [2.47, 2.47, 2.47, 2.47, 2.47, 2.58, 2.58, 2.63]  # g/cm3

    # dipsersion parameters
    def f(): return np.logspace(np.log10(0.2), np.log10(3.5), 85)
    curves = [('R', 'U', 0, f()),
              ('R', 'U', 1, f()),
              ('R', 'C', 0, f()),
              ('R', 'C', 1, f()),
              ('L', 'U', 0, f()),
              ('L', 'U', 1, f()),
              ('L', 'C', 0, f()),
              ('L', 'C', 1, f())]

    # reorganize inputs for dispersion_1
    Waves, Types, Modes, Freqs = zip(*curves)

    # compute dispersion curves
    with Timer('dispersion'):
        out = list(dispersion_1(ztop, vp, vs, rh, Waves, Types, Modes, Freqs))

    # display results
    ax = plt.gca()
    for w, t, m, fs, us in out:
        ax.loglog(1. / fs, us, '+-', label="%s%s%d" % (w, t, m))
    ax.set_xlabel('period (s)')
    ax.set_ylabel('velocity (km/s)')
    ax.grid(True, which="major")
    ax.grid(True, which="minor")
    logtick(ax, "xy")
    plt.legend()
    plt.show()
