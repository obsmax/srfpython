#!/usr/bin/python2.7
from subprocess import Popen, PIPE
import numpy as np
import sys, glob, os, time, signal

_pathfile = os.path.realpath(__file__) #.../srfpyhon/src/srfdis17.py
_file     = _pathfile.split('/')[-1]
srfpre96_exe = _pathfile.replace('/src/', '/bin/').replace(_file, 'max_srfpre96') #.../srfpyhon/bin/max_srfpre96
srfdis96_exe = _pathfile.replace('/src/', '/bin/').replace(_file, 'max_srfdis96') #.../srfpyhon/bin/max_srfdis96
if not os.path.exists(srfpre96_exe) or not os.path.exists(srfdis96_exe):
    raise Exception('could not find %s and/or %s' % (srfpre96_exe, srfdis96_exe))

###################################### TOOLS
class Timeout():
    class Timeout(Exception): pass
    def __init__(self, sec):
        self.sec = sec
    def __enter__(self):
        signal.signal(signal.SIGALRM, self.raise_timeout)
        signal.alarm(self.sec)
    def __exit__(self, *args):
        signal.alarm(0)
    def raise_timeout(self, *args):
        raise Timeout.Timeout()
#_____________________________________
def firstfalse(I):
    if not I[0] : return 0
    II = I[1:] != I[:-1]
    r = np.arange(len(I) - 1)
    return r[II][0] + 1
#_____________________________________
def munique(*Xs):
    assert np.all(len(X) == len(Xs[0]) for X in Xs)
    L = []
    for tup in zip(*Xs):
        if tup not in L: 
            L.append(tup)
    return tuple([np.array(w) for w in zip(*L)])
#_____________________________________
def freqspace(freqmin, freqmax, nfreq, scale="flin"):
    if "lin" in scale.lower():
        return np.linspace(freqmin, freqmax, nfreq)
    elif "log" in scale.lower():
        return np.logspace(np.log10(freqmin), np.logspace(freqmax), nfreq)
    else: raise ValueError('%s not understood' % scale)
###################################### MODIFIED HERRMANN'S CODES
def prep_srfpre96_1(h = 0.005, dcl = 0.005, dcr = 0.005):
    """prepare input for modified srfpre96 (max_srfpre96)
    dcl, dcr are root search increment for loev and rayleigh waves respectively
    h, dcl,dcr : see srfpre96
    """
    return "%f %f %f" % (h, dcl, dcr)
#_____________________________________
def prep_srfpre96_2(z, vp, vs, rh):
    """prepare input for modified srfpre96 (max_srfpre96)
       z   = depth in km, z[0] must be 0
       vp  = vp in km/s
       vs  = vs in km/s
       rh  = density in g/cm3
    """

    if z[0]: raise CPiSDomainError('Z0 must be 0')#assert z[0] == 0.
    n = len(z)

    if not (n == len(vp) == len(vs) == len(rh)): 
        raise CPiSDomainError('Z VP, VS, RH must be the same length')#    assert n == len(vp) == len(vs) == len(rh)
    z = np.asarray(z, float)
    if (np.isinf(z) | np.isnan(z)).any():
        raise CPiSDomainError('Z must be growing, layers must be at least 0.001km thick, got %s' % str(z))
    if (z[1:] - z[:-1] < 0.001).any():#    assert np.all(z[1:] - z[:-1] >= 0.001)
        raise CPiSDomainError('Z must be growing, layers must be at least 0.001km thick, got %s' % str(z))


    vs = np.asarray(vs, float)
    if (np.isinf(vs) | np.isnan(vs)).any(): raise CPiSDomainError('vs value error %s' % str(vs))
    if not (vs > 0.08).all():               raise CPiSDomainError('vs domain error %s' % str(vs))

    vp = np.asarray(vp, float)
    if (np.isinf(vp) | np.isnan(vp)).any():      raise CPiSDomainError('vp value error %s' % str(vp))
    if not ((vp > vs * (4. / 3.) ** 0.5)).all(): raise CPiSDomainError('vp over vs domain error %s' % str(vp / vs))

    rh = np.asarray(rh, float)
    if (np.isinf(rh) | np.isnan(rh)).any():      raise CPiSDomainError('rh value error %s' % str(rh))
    if not  np.all((rh > 1.)):                   raise CPiSDomainError('density domain error : %s' % str(rh))

    out="""MODEL.01
model-title
ISOTROPIC
KGS
FLAT EARTH
1-D
CONSTANT VELOCITY
LINE08
LINE09
LINE10
LINE11
      H(KM)   VP(KM/S)   VS(KM/S) RHO(GM/CC)     QP         QS       ETAP       ETAS      FREFP      FREFS    """
    fmt = "\n %f %f %f %f 0.00 0.00 0.00 0.00 1.00 1.00"

    for i in xrange(n - 1):
        out += fmt % (z[i+1] - z[i], vp[i], vs[i], rh[i])
    out += fmt % (0., vp[n-1], vs[n-1], rh[n-1])
    return out
#_____________________________________
def prep_srfpre96_3(waves,types,modes,freqs):
    """prepare input for modified srfpre96 (max_srfpre96)"""
    nlc = nlu = nrc = nru = 0
    fmt="\nSURF96 %1s %1s X   %d  %f 1. 1."
    out = ""
    for w,t,m,f in zip(waves,types,modes,freqs):
        out += fmt % (w,t,m,1./f)
        if w == "R":
            if    t == "C": nrc = max([nrc, m + 1])
            elif  t == "U": nru = max([nru, m + 1])
        elif w == "L":
            if    t == "C": nlc = max([nlc, m + 1])
            elif  t == "U": nlu = max([nlu, m + 1])

    out = "%d %d %d %d" %(nlc,nlu,nrc,nru) + out
    return out
#_____________________________________
def groupbywtm(waves, types, modes, freqs, values, dvalues = None, keepnans = True):
    """group outputs from srfdis17 by wave, type, modes

    waves types modes freqs are the same length
    groupbywtm demultiplexes these arrays to produce, for each mode
    w (wave key, scalar), t (type key, scalar), m (mode number, scalar), freqs (array), values (array), dvalues (array)
    
    keepnans : by default the nan are removed from the arrays (this occurs especially for overtones below the cut-off period)

    """
    freqs = np.array(freqs, float)
    types = np.array(list(types), "|S1")
    waves = np.array(list(waves), "|S1")
    modes = np.array(modes,  int)
    assert len(waves) == len(types) == len(modes) == len(freqs)
    if dvalues is not None:
        dvalues = np.array(dvalues, float)
        assert len(waves) == len(dvalues)

    w, t, m = munique(waves, types, modes)
    I = np.lexsort((m, t, w))

    for w, t, m in zip(w[I], t[I], m[I]):
        J = (waves == w) & (types == t) & (modes == m)
        
        if keepnans:   K = np.ones(len(values[J]), bool)
        else:          K = ~np.isnan(values[J])
        
        L = np.argsort(freqs[J][K])

        if dvalues is None:
            yield w, t, m, freqs[J][K][L], values[J][K][L]
        else:
            yield  w, t, m, freqs[J][K][L], values[J][K][L], dvalues[J][K][L]
#_____________________________________
def igroupbywtm(Waves, Types, Modes, Freqs):
    """
    make the opposite of groupbywtm, prepare input for srfdis17
    e.g. 
    f = freqspace(0.1, 20, 5, 'plog')
    Waves                = ['R', 'R', 'R', 'R']
    Types                = ['U', 'U', 'C', 'C']
    Modes                = [ 0,   1,   0,   1 ]
    Freqs                = [ f,   f,   f,   f ]
    
    becomes 
    array(['R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R','R', 'R', 'R', 'R', 'R', 'R', 'R'], dtype='|S1')
    array(['U', 'U', 'U', 'U', 'U', 'U', 'U', 'U', 'U', 'U', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'], dtype='|S1')
    array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
    array([0.1, 0.37606031, 1.41421356, 5.3182959, 20., 0.1, 0.37606031, 1.41421356, 5.3182959, 20., 0.1, 0.37606031, 1.41421356, 5.3182959, 20., 0.1, 0.37606031, 1.41421356, 5.3182959, 20])
    """
    Waves, Types, Modes = [np.array(w, object) for w in Waves, Types, Modes]
    waves = Waves.repeat([len(ff) for ff in Freqs])
    types = Types.repeat([len(ff) for ff in Freqs])
    modes = Modes.repeat([len(ff) for ff in Freqs])
    waves = np.array(waves, '|S1')
    types = np.array(types, '|S1')
    modes = np.array(modes, int)
    freqs = np.concatenate(Freqs)
    return waves, types, modes, freqs 
#_____________________________________
class CPiSError(Exception): pass
class CPiSDomainError(CPiSError): pass
#_____________________________________
def readsrfdis96(stdout, waves, types, modes, freqs):
    "converts output from max_srfdis96, to be optimized (takes about 0.5 * the times used by max_srfpre96+max_srfdis96)"

    #try:     A = np.genfromtxt(StringIO(unicode(stdout)), dtype = float)
    #except:  raise CPiSError('could not understand stdout \n%s' % stdout)
    A = np.asarray( (" ".join(stdout.strip().rstrip('\n').split('\n'))).split(), float)
    A = A.reshape((len(A) / 6, 6))

    W = A[:, 0]
    M = A[:, 1]
    I = A[:, 3] == 0.
    nI = ~I
    n = A.shape[0]
    T, C, U = np.zeros(n, float), np.zeros(n, float), np.zeros(n, float) * np.nan
    T[I]=A[I,2]
    T[nI]=A[nI,2:4].mean(axis = 1)
    C[I]=A[I,4]
    C[nI]=A[nI,4:6].mean(axis = 1)

    U[nI] = C[nI] / (1. - (np.log(A[nI,5]) - np.log(A[nI,4])) / (np.log(A[nI,2]) - np.log(A[nI,3])))
    
    L, R = (W == 1), (W == 2)
    Ms = {}
    for m in np.unique(modes): Ms['%d' % m] = (M == m)


    values = np.zeros(len(waves), float) * np.nan
    for n, (w, t, m, f) in enumerate(zip(waves, types, modes, freqs)):
        #S = eval(w) & (M == m)
        if   w == "R":S = R & Ms['%d' % m]
        elif w == "L":S = L & Ms['%d' % m]

        if not S.any(): continue#yield w, t, int(m), f, np.nan
        TS = T[S]
        iS = np.abs(TS - 1./f).argmin()
        per = TS[iS]        
        if abs(per - 1./f) * f > 0.01: continue#yield w, t, int(m), f, np.nan

        if   t == 'C': val = C[S][iS] 
        elif t == "U": val = U[S][iS]
        else: raise ValueError('')
        #if val <= 0: #NON, je met une penalite dans la fonction cout
        #    raise CPiSError('encountered negative dispersion velocity')

        values[n] = val

    return values



#_____________________________________
def srfdis17(ztop, vp, vs, rh, \
    waves, types, modes, freqs,
    h = 0.005, dcl = 0.005, dcr = 0.005):

    """srfdis17 : compute dispersion curves from a depth model for desired wave (R or L) types (C or U for phase or group) and frequencies (Hz)
                  based on a modified version of the codes from Herrmann and Ammon 2002

    input: 
        -> depth model
        ztop  : list or array, top layer depth in km, ztop[0] must be 0 !!!
        vp    : list or array, km/s
        vs    : list or array, km/s
        rh    : list or array, density in g.cm^-3
                note that len(ztop) must be len(vp) = len(vs) = len(rh)

        -> required dispersion points
        waves : list or array, like ['L', 'L', 'L', 'R', 'R', 'R', 'R', 'R'] (L = Love, R = Rayleigh)
        types : list or array, like ['U', 'U', 'U', 'C', 'C', 'C', 'C', 'C'] (U = groupe, C = phase)
        modes : list or array, like [ 0,   0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ] (mode number, 0 = fundamental)
        freqs : list or array, like [ 1.,  2.,  3.,  1.,  2.,  3.,  4.,  5.] (frequency in Hz)
                note that len(waves) must be len(types) = len(modes) = len(freqs)

        -> Herrmann's parameters, see CPS documentation
        h = period increment for phase to group conversion (0.005 is reasonale)
        dcl, dcr = phase velocity increment for love and rayleigh root searches respectively, see Herrmann doc

    output:
        values : list or array, dispersion values for each wave, type, mode, possibly nan (cut-off period)
                 note : use groupbywtm to group outputs by wave, type, and mode
    see also :
        srfdis17_1
        groupbywtm
        igroupbywtm
    """

    #-------------- 
    instr1 = prep_srfpre96_1(h = h, dcl = dcl, dcr = dcr)
    instr2 = prep_srfpre96_2(ztop, vp, vs, rh)
    instr3 = prep_srfpre96_3(waves, types, modes, freqs)
    pstdin = "\n".join([instr1, instr2, instr3])

    #-------------- 
    try:
        with Timeout(5):
            p = Popen(srfpre96_exe, stdout=PIPE, stdin=PIPE, stderr=PIPE, shell=False, preexec_fn=os.setsid)#, stderr = stderrfid)
            p.stdin.write(pstdin)
            pstdout, _ = p.communicate()
            if p.returncode : raise CPiSError('srfpre96_exe failed')

            q = Popen(srfdis96_exe, stdout=PIPE, stdin=PIPE, stderr = PIPE, shell = False,  preexec_fn=os.setsid)#, stderr = stderrfid)
            q.stdin.write(pstdout)
            qstdout, _ = q.communicate()
            if q.returncode : raise CPiSError('srfdis96_exe failed')    
    except Timeout.Timeout:
        print "*123*", ztop, vp, vs, rh
        #raise Exception('srfdis96 timed out ')
        os.killpg(os.getpgid(p.pid), signal.SIGKILL)
        os.killpq(os.getpgid(q.pid), signal.SIGKILL)
        raise CPiSError('timeout')
    except: raise 
    finally:
        p.stdin.close()
        q.stdin.close()

    #-------------- 
    values = readsrfdis96(qstdout, waves, types, modes, freqs)
    return values

#_____________________________________
def srfdis17_1(ztop, vp, vs, rh, \
    Waves, Types, Modes, Freqs,
    h = 0.005, dcl = 0.005, dcr = 0.005, keepnans = False):

    """same as srfdis17 with slightely more convenient input and output (no need to repeat wave, type and mode)

    Waves is like ['L', 'L', 'R']
    Types is like ['C', 'C', 'U']
    Modes is like [ 0,   1,   0 ]
    Freqs is like [fLC0, fLC1, fRU0] where f??? are frequency numpy arrays or lists

    output is like
        'L', 'C', 0, fLC0, LC0(fLC0)
        'L', 'C', 1, fLC1, LC1(fLC1)
        'R', 'U', 0, fRU0, RU0(fRU0)

    """
    waves, types, modes, freqs = igroupbywtm(Waves, Types, Modes, Freqs)
    values = srfdis17(ztop, vp, vs, rh, \
                        waves, types, modes, freqs,
                        h = h, dcl = dcl, dcr = dcr)

    for w, t, m, F, V in groupbywtm(waves, types, modes, freqs, values, keepnans = keepnans):
        yield w, t, m, F, V
#_____________________________________
if __name__ == "__main__":
    """ DEMO """
    import matplotlib.pyplot as plt

    ###depth model
    ztop = [0.00, 0.25, 0.45, 0.65, 0.85, 1.05, 1.53, 1.80] #km, top layer depth
    vp   = [1.85, 2.36, 2.63, 3.15, 3.71, 4.54, 5.48, 5.80] #km/s
    vs   = [0.86, 1.10, 1.24, 1.47, 1.73, 2.13, 3.13, 3.31] #km/s
    rh   = [2.47, 2.47, 2.47, 2.47, 2.47, 2.58, 2.58, 2.63] #g/cm3

    ###dipsersion parameters
    f = np.logspace(np.log10(0.2), np.log10(10.), 20)
    Waves = ['R', 'R', 'R', 'R']
    Types = ['U', 'U', 'U', 'U']
    Modes = [ 0 ,  1,   2 ,  3 ]
    Freqs = [ f ,  f,   f,   f ]

    start = time.time()
    out = list(srfdis17_1(ztop, vp, vs, rh, Waves, Types, Modes, Freqs))
    print "elapsed time %fs" % (time.time() - start)
    for w, t, m, fs, vs in out:
        plt.gca().loglog(1. / fs, vs, '+-')
    plt.gcf().show()
    raw_input('pause')

