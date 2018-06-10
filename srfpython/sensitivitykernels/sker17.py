#!/usr/bin/env python

import numpy as np
from srfpython.Herrmann.Herrmann import dispersion, dispersion_1, Timer, groupbywtm, igroupbywtm, CPiSDomainError
from srfpython.utils import minmax
from srfpython.standalone.cmaps import cccfcmap3, tomocmap1
from srfpython.standalone.multipro8 import MapSync, Job
from srfpython.standalone.display import logtick

"""
srfker17, Maximilien Lehujeur, 01/11/2017
module to compute finite difference surface wave sensitivity kernels 
see documentation in function sker17
use __main__ for demo

see also Herrmann.dispersion.dispersion
"""


# _____________________________________
def lognofail(x):

    def ilognofail(x):
        if np.isnan(x):  return x
        elif x < 0.:     return np.nan
        elif x == 0.:    return -np.inf
        return np.log(x)

    if hasattr(x, "__iter__"):
        return np.asarray(map(ilognofail, x), float)
    else:
        return ilognofail(x)


# _____________________________________
def sker17(ztop, vp, vs, rh, \
    waves, types, modes, freqs,
    dz=0.001, dlogvs=0.01, dlogpr=0.01, dlogrh=0.01, norm=True,
    h = 0.005, dcl = 0.005, dcr = 0.005):
    """sker17 : compute finite difference sensitivity kernels for surface waves dispersion curves 
    input: 
        -> depth model
        ztop, vp, vs, rh  : lists or arrays, see dispersion

        -> required dispersion points
        waves, types, modes, freqs : lists or arrays, see dispersion

        -> sensitivity kernel computation
        dz = depth increment in km
        dlogvs = increment to apply to the logarithm of vs
        dlogpr = increment to apply to the logarithm of vp/vs
        dlogrh = increment to apply to the logarithm of rho
        norm = if True, I divide the sensitivity values by the thickness of each layer
                => this corrects for the difference of sensitivity due to the variable thicknesss

        -> Herrmann's parameters, see CPS documentation
        h, dcl, dcr = passed to dispersion

    output:
        -> yields a tuple (w, t, m, F, DLOGVADZ, DLOGVADLOGVS, DLOGVADLOGPR, DLOGVADLOGRH) for each wave, type and mode
        w      = string, wave letter (L = Love or R = Rayleigh)
        t      = string, type letter (C = phase or U = group)
        m      = int, mode number (0= fundamental)
        F      = array, 1D, frequency array in Hz
        DLOGVADZ  = array, 2D, [normed] sensitivity kernel relative to top depth of each layer (lines) and frequency (columns)
        DLOGVADLOGVS  = array, 2D, [normed] sensitivity kernel relative to Pwave velocity of each layer (lines) and frequency (columns)
        DLOGVADLOGPR  = array, 2D, [normed] sensitivity kernel relative to Swave velocity of each layer (lines) and frequency (columns)
        DLOGVADLOGRH  = array, 2D, [normed] sensitivity kernel relative to density of each layer (lines) and frequency (columns)                
                 note that these arrays might contain nans
    see also :
        sker17_1
        dispersion
    """

    waves, types, modes, freqs = [np.asarray(_) for _ in waves, types, modes, freqs]
    nlayer = len(ztop)
    H = np.array(ztop) # NOT ASARRAY
    H[:-1], H[-1] = H[1:] - H[:-1], np.inf #layer thickness in km

    model0 = np.concatenate((ztop, np.log(vs), np.log(vp/vs), np.log(rh)))
    dmodel = np.concatenate((dz * np.ones_like(ztop),
                             dlogvs * np.ones_like(vs),
                             dlogpr * np.ones_like(vs),
                             dlogrh * np.ones_like(rh)))

    logvalues0 = lognofail(dispersion(ztop, vp, vs, rh,
                       waves, types, modes, freqs,
                       h = h, dcl = dcl, dcr = dcr))

    IZ = np.arange(nlayer)
    IVS = np.arange(nlayer, 2*nlayer)
    IPR = np.arange(2*nlayer, 3*nlayer)
    IRH = np.arange(3*nlayer, 4*nlayer)        
    DVADP = np.zeros((4 * nlayer, len(waves)), float) * np.nan

    # ----
    # parallel
    # ----
    def fun(i, modeli):
        ztopi, logvsi, logpri, logrhi = \
            modeli[IZ], modeli[IVS], modeli[IPR], modeli[IRH]
        n = len(ztopi)
        ilayer = i % n
        if ilayer == n-1:
            Hi = 1.e50 # thickness of the half-space
        else:
            Hi = ztopi[ilayer + 1] - ztopi[ilayer]

        try:
            logvaluesi = lognofail(dispersion(
                ztopi, np.exp(logvsi + logpri),
                np.exp(logvsi), np.exp(logrhi),
                waves, types, modes, freqs,
                h=h, dcl=dcl, dcr=dcr))
        except CPiSDomainError as err:
            print ("error during gradient computation %s" % str(err))
            return i, None
        except:
            raise
        if norm:
            # sensitivity corrected from the layer thicknesses
            DVAVPi = (logvaluesi - logvalues0) / (modeli[i] - model0[i]) / Hi
        else:
            # absolute sensitivity regardless the thickness differences
            DVAVPi = (logvaluesi - logvalues0) / (modeli[i] - model0[i])

        return i, DVAVPi

    # ----
    def gen():
        for i in xrange(1, 4 * len(ztop)):
            modeli = model0.copy()
            modeli[i] += dmodel[i]
            yield Job(i, modeli)

    # ----
    with MapSync(fun, gen()) as ma:
        for _, (i, DVAVPi), _, _ in ma:
            if DVAVPi is None: continue
            DVADP[i, :] = DVAVPi

    for w, t, m, F, Iwtm in groupbywtm(waves, types, modes, freqs, np.arange(len(waves))):
        DLOGVADZ = DVADP[IZ, :][:, Iwtm]
        DLOGVADLOGPR = DVADP[IPR, :][:, Iwtm]
        DLOGVADLOGVS = DVADP[IVS, :][:, Iwtm]
        DLOGVADLOGRH = DVADP[IRH, :][:, Iwtm]
        DLOGVADZ, DLOGVADLOGVS, DLOGVADLOGPR, DLOGVADLOGRH = \
            [np.ma.masked_where(np.isnan(_), _) for _ in
             [DLOGVADZ, DLOGVADLOGVS, DLOGVADLOGPR, DLOGVADLOGRH]]
        
        yield w, t, m, F, DLOGVADZ, DLOGVADLOGVS, DLOGVADLOGPR, DLOGVADLOGRH


# _____________________________________
def sker17_1(ztop, vp, vs, rh,
    Waves, Types, Modes, Freqs, **kwargs):
    """sker17_1 : same as sker17 with slightely more convenient input (no need to repeat wave, type and mode)

    Waves is like ['L', 'L', 'R']
    Types is like ['C', 'C', 'U']
    Modes is like [ 0,   1,   0 ]
    Freqs is like [fLC0, fLC1, fRU0] where f??? are frequency numpy arrays or lists
    
    see sker17 for detailed input and output arguments
    """
    waves, types, modes, freqs = igroupbywtm(Waves, Types, Modes, Freqs)
    for tup in sker17(ztop, vp, vs, rh,
            waves, types, modes, freqs, **kwargs):
        yield tup


# -----------------------------
if __name__ == "__main__":
    import sys

    help = '''sker17
    -m96          depthmodel to read 
    -RU0          rayleigh, group, mode 0 : expects 4 frequency arguments : fstart, fend, nfreq, fscale
    -RU1          rayleigh, group, mode 1 : expects 4 frequency arguments : fstart, fend, nfreq, fscale
    -RC0          rayleigh, phase, mode 0 : expects 4 frequency arguments : fstart, fend, nfreq, fscale          
    -LC0          love,     phase, mode 0 : expects 4 frequency arguments : fstart, fend, nfreq, fscale
    ...
    -png          save figures as pngfiles (overwrite if exists)
                  sker17_depthdisp.png
                  sker17_RU0_fstart_fend_nfreq_fscale.png
                  ...
    '''
    if len(sys.argv) == 1:
        print help
        sys.exit()

    from srfpython.standalone.display import plt
    from srfpython.utils import readargv
    from srfpython.Herrmann.Herrmann import dispersion_1
    from srfpython.depthdisp.dispcurves import freqspace
    from srfpython.depthdisp.depthmodels import depthmodel_from_mod96
    import numpy as np

    argv = readargv()
    # -----------------------------------:
    dm = depthmodel_from_mod96(argv['m96'][0])
    ztop = dm.vp.ztop()
    vp = dm.vp.values
    vs = dm.vs.values
    rh = dm.rh.values

    # -----------------------------------
    Waves, Types, Modes, Freqs = [], [], [], []
    for k in argv.keys():
        if k[0].upper() in "RL" and k[1].upper() in "UC" and k[2] in "0123456789":
            fstart, fend, nfreq, fspace = argv[k]
            freq = freqspace(float(fstart), float(fend), int(nfreq), fspace)
            Waves.append(k[0])
            Types.append(k[1])
            Modes.append(int(k[2:]))
            Freqs.append(freq)

    # -----------------------------------
    # ##compute dispersion curves
    fig1 = plt.figure(figsize=(8,4))
    fig1.subplots_adjust(wspace=0.3)
    with Timer('dispersion'):
        out = list(dispersion_1(ztop, vp, vs, rh, Waves, Types, Modes, Freqs))
    ax1 = fig1.add_subplot(121)
    dm.show(ax1)
    ax1.grid(True, linestyle=":", color="k")
    plt.legend()
    ax2 = fig1.add_subplot(122)
    for w, t, m, fs, us in out:
        ax2.loglog(1. / fs, us, '+-', label="%s%s%d" % (w, t, m))
    ax2.set_xlabel('period (s)')
    ax2.set_ylabel('velocity (km/s)')
    ax2.grid(True, which="major")
    ax2.grid(True, which="minor")
    logtick(ax2, "xy")
    plt.legend()
    if "png" in argv.keys():
        fout = "sker17_depthdisp.png"
        print fout
        fig1.savefig(fout, dpi=300)
    else:
        fig1.show()


    # ##sensitivity kernels
    norm = True
    fig = plt.figure()
    fig.subplots_adjust(wspace=0.1, hspace=0.2)
    if "png" not in argv.keys():
        fig.show()

    for w, t, m, F, DLOGVADZ, DLOGVADLOGVS, DLOGVADLOGPR, DLOGVADLOGRH in \
            sker17_1(ztop, vp, vs, rh,
                     Waves, Types, Modes, Freqs,
                     dz=0.001, dlogvs=.01, dlogpr=.01, dlogrh=.01, norm=True,
                     h=0.005, dcl=0.005, dcr=0.005):

        fig.clf()
        # ------
        _depth_ = np.concatenate((ztop, [1.1 * ztop[-1]]))
        _F_ = np.concatenate(([F[0] * 0.95], F * 1.05))
        fig.suptitle('%s%s%d' % (w, t, m))

        # ------
        vmax = abs(DLOGVADLOGVS).max()#np.max([abs(DLOGVADZ).max(), abs(DLOGVADLOGVS).max(), abs(DLOGVADLOGPR).max(), abs(DLOGVADLOGRH).max()])
        vmax = np.min([vmax, 10.])
        vmin, vmax, cmap = -vmax, vmax, tomocmap1(W=.25) #cccfcmap3() #plt.cm.RdBu
        ax1 = fig.add_subplot(221)
        plt.pcolormesh(1. / _F_, _depth_, DLOGVADZ,
                       vmin=vmin, vmax=vmax, cmap=cmap)

        ax2 = fig.add_subplot(222, sharex=ax1, sharey=ax1)
        plt.pcolormesh(1. / _F_, _depth_, DLOGVADLOGVS,
                       vmin=vmin, vmax=vmax, cmap=cmap)

        ax3 = fig.add_subplot(223, sharex=ax1, sharey=ax1)
        plt.pcolormesh(1. / _F_, _depth_, DLOGVADLOGPR,
                       vmin=vmin, vmax=vmax, cmap=cmap)

        ax4 = fig.add_subplot(224, sharex=ax1, sharey=ax1)
        plt.pcolormesh(1. / _F_, _depth_, DLOGVADLOGRH,
                       vmin=vmin, vmax=vmax, cmap=cmap)

        cax = fig.add_axes((.91, .3, .01, .4))
        cb = plt.cm.ScalarMappable(norm=None, cmap=cmap)
        cb.set_array([vmin, vmax])
        fig.colorbar(cb, cax=cax)

        #  ------
        for ax, p in zip([ax1, ax2, ax3, ax4], ["Z^{top}_i", "ln Vs_i", "ln (Vp/Vs)_i", r"ln \rho _i"]):
            ax.set_xlim(minmax(1. / _F_))
            ax.set_ylim(minmax(_depth_))
            ax.set_xscale('log')
            if ax in [ax1, ax3]:
                ax.set_ylabel('depth (km)')
            if ax in [ax3, ax4]:
                ax.set_xlabel('period (s)')

            ax.set_title(r'$ \frac{1}{H_i} \, \frac{d ln%s_j}{d %s} $' % (t, p))

            if not ax.yaxis_inverted(): ax.invert_yaxis()
            logtick(ax, "x")

        # ------
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.setp(ax4.get_yticklabels(), visible=False)
        fig.canvas.draw()

        if "png" in argv.keys():
            k = "%s%s%d" % (w, t, m)
            fout = 'sker17_%s_%s_%s_%s_%s.png' % (k, argv[k][0], argv[k][1], argv[k][2], argv[k][3])
            print fout
            fig.savefig(fout, dpi=300)
        else:
            raw_input('pause : press enter to plot the next wave type and mode')
    # --------------------
    if "png" not in argv.keys():
        raw_input('bye')

