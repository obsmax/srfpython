#!/usr/bin/env python
from typing import List
from builtins import input
from dataclasses import dataclass
import sys
import numpy as np

from srfpython.depthdisp.dispcurves import freqspace
from srfpython.depthdisp.depthmodels import depthmodel, depthmodel_from_mod96
from srfpython.Herrmann.Herrmann import Timer, groupbywtm, igroupbywtm, CPiSDomainError
from srfpython.Herrmann.Herrmann import HerrmannCallerBasis, HerrmannCallerFromGroupedLists
from srfpython.Herrmann.Herrmann import Curve
from srfpython.utils import minmax, readargv
from srfpython.standalone.cmaps import cccfcmap3, tomocmap1
from srfpython.standalone.multipro8 import MapSync, Job
from srfpython.standalone.display import logtick, textonly
from srfpython.standalone.display import plt

"""
srfker17, Maximilien Lehujeur, 01/11/2017
module to compute finite difference surface wave sensitivity kernels 
see documentation in function sker17
use __main__ for demo

see also Herrmann.dispersion.dispersion
"""


def lognofail(x):

    def ilognofail(x):
        if np.isnan(x):  return x
        elif x < 0.:     return np.nan
        elif x == 0.:    return -np.inf
        return np.log(x)

    if hasattr(x, "__iter__"):
        return np.asarray(list(map(ilognofail, x)), float)
    else:
        return ilognofail(x)


@dataclass
class SensitivityKernel:
    parameter_name: str = None
    norm: bool = None
    relative: bool = None
    curve: Curve = None
    dm: depthmodel = None
    kernel: np.ndarray = None
    all_curves: List[Curve] = None  # keep also all other curves for display

    def expression(self):
        param_latex = {
            "ztop": r"Z^{top}",
            "vp": r"Vp",
            "vs": r"Vs",
            "rh": r"\rho",
            }[self.parameter_name]
        if self.relative:
            if self.norm:
                txt = r'$ \frac{H}{H_i} \, \frac{d ln%s_j}{d ln%s_i} $' % (self.curve.type, param_latex)
            else:
                txt = r'$ \frac{d ln%s_j}{d ln%s_i} $' % (self.curve.type, param_latex)
        else:
            if self.norm:
                txt = r'$ \frac{H}{H_i} \, \frac{d %s_j}{d %s_i} $' % (self.curve.type, param_latex)
            else:
                txt = r'$ \frac{d %s_j}{d %s_i} $' % (self.curve.type, param_latex)
        return txt

    @property
    def z_edges(self):
        return np.hstack((self.dm.vp.z, [1.1 * self.dm.vp.z[-1]])) # km

    @ property
    def z_mid(self):
        return self.dm.vp.z_mid  # km

    @property
    def f_edges(self):
        return np.hstack((
            self.curve.freqs[0] * 0.95,
            np.sqrt(self.curve.freqs[1:] * self.curve.freqs[:-1]),
            self.curve.freqs[-1] * 1.05))  # Hz

    def show(self, units: str='m/s/kg.m-3',
             vmin: float=None, vmax: float=None,
             cmap=None):

        if units == 'm/s/kg.m-3':
            kf = 1.0
            kx = 1000.
            kv = 1000.
        else:
            raise NotImplemented

        # ==== display
        fig1 = plt.figure(figsize=(8, 8))
        fig1.subplots_adjust(wspace=0.3)
        gs = fig1.add_gridspec(3, 3)

        # depth model
        ax1 = fig1.add_subplot(gs[1:, 0])
        self.dm.show(ax1, ".-", units=units, )
        ax1.grid(True, linestyle=":", color="k")
        plt.legend()

        # disp curve
        ax2 = fig1.add_subplot(
            gs[0, 1:],
            ylabel={'m/s/kg.m-3': 'Velocity (m/s)'}[units])

        for curve in self.all_curves:
            ax2.plot(
                curve.freqs * kf,
                curve.values * kv,
                '+-',
                label=f"{curve.wave}{curve.type}{curve.mode}")
        #
        # ax2.plot(
        #     self.curve.freqs * kf,
        #     self.curve.values * kv,
        #     '+-',
        #     label=f"{self.curve.wave}{self.curve.type}{self.curve.mode}")

        ax2.grid(True, which="major")
        ax2.grid(True, which="minor")

        if False:
            ax2.set_xscale('log')
            ax2.set_yscale('log')
            logtick(ax2, "xy")

        plt.legend()

        ax3 = fig1.add_subplot(gs[1:, 1:], sharex=ax2, sharey=ax1)
        ax3.set_title(f"{self.curve.wave}{self.curve.type}{self.curve.mode}")

        if vmax is None:
            vmax = abs(self.kernel).max()
        if vmin is None:
            vmin = -vmax

        coll = ax3.pcolormesh(
            self.f_edges * kf,
            self.z_edges * kx,
            self.kernel,
            vmin=vmin, vmax=vmax, cmap=cmap)

        if False:
            if skernel.max() - skernel.min():
                levels = np.logspace(-1., 2, 10)
                levels = np.hstack((-levels[::-1], [0], levels))
                plt.contour(freqs, z_mid, skernel,
                            levels=levels,
                            colors="k")

        cax = fig1.add_axes((.91, .2, .01, .2))
        plt.colorbar(coll, cax=cax)

        ax3.set_xlabel({'m/s/kg.m-3': 'Frequency (Hz)'}[units])

        ax3.set_xlim(min(self.f_edges * kf), max(self.f_edges * kf))
        ax3.set_ylim(min(self.z_edges * kx), max(self.z_edges * kx))
        if False:
            ax3.set_xscale('log')
            logtick(ax3, "x")

        textonly(ax3, self.expression(), loc=3, fontsize=16)

        if not ax3.yaxis_inverted():
            ax3.invert_yaxis()

        return fig1


def sker17(
        # ztop: np.ndarray, vp: np.ndarray, vs: np.ndarray, rh: np.ndarray,
        dm: depthmodel,
        # waves: np.ndarray, types: np.ndarray, modes: np.ndarray, freqs: np.ndarray,
        curves: List[Curve],
        dlogz: float=0.001, dlogvp: float=0.01, dlogvs: float=0.01, dlogrh: float=0.01,
        norm: bool=True,
        relative: bool=True,
        h: float=0.005, ddc: float=0.005,
        ):

    """sker17 : compute finite difference sensitivity kernels for surface waves dispersion curves 
    input: 
        -> depth model
        ztop, vp, vs, rh  : lists or arrays, see dispersion

        -> required dispersion points
        waves, types, modes, freqs : lists or arrays, see dispersion

        -> sensitivity kernel computation
        dlogz  = relative z increment, dimless
        dlogvp = relative vp increment, dimless
        dlogvs = relative vs increment, dimless
        dlogrh = relative rh increment, dimless
        norm:
            if False: no layer thickness correction applied, the output kernels are
                    dlnc)fj / dztop)zi
                    dlnc)fj / dlnvp)zi
                    dlnc)fj / dlnvs)zi
                    dlnc)fj / dlnrh)zi
            if True: The output kernels are multiplied by H / Hi
                so that heterogeneous layer thickness are taken into account and the dimension of the kernels is preserved

        -> Herrmann's parameters, see surf96 documentation (CPS)
        h, ddc = passed to dispersion

    output:
        -> yields a tuple (w, t, m, F, DLOGVADZ, DLOGVADLOGVS, DLOGVADLOGPR, DLOGVADLOGRH) for each wave, type and mode
        w      = string, wave letter (L = Love or R = Rayleigh)
        t      = string, type letter (C = phase or U = group)
        m      = int, mode number (0= fundamental)
        F      = array, 1D, frequency array in Hz
        DLOGVADZ  = array, 2D, [normalized] sensitivity kernel relative to top depth of each layer (lines) and frequency (columns)
        DLOGVADLOGVP  = array, 2D, [normalized] sensitivity kernel relative to Pwave velocity of each layer (lines) and frequency (columns)
        DLOGVADLOGVS  = array, 2D, [normalized] sensitivity kernel relative to Swave velocity of each layer (lines) and frequency (columns)
        DLOGVADLOGRH  = array, 2D, [normalized] sensitivity kernel relative to density of each layer (lines) and frequency (columns)                
                 note that these arrays might contain nans
    see also :
        sker17_1
        dispersion
    """

    Waves = [curve.wave for curve in curves]
    Types = [curve.type for curve in curves]
    Modes = [curve.mode for curve in curves]
    Freqs = [curve.freqs for curve in curves]

    waves, types, modes, freqs = igroupbywtm(Waves=Waves, Types=Types, Modes=Modes, Freqs=Freqs)
    herrmanncaller = HerrmannCallerBasis(waves, types, modes, freqs, h=h, ddc=ddc)

    nlayer = dm.vp.z.size
    H = dm.vp.thickness()   # layer thickness in km, ends by np.inf

    # reference model for differentiation
    model0 = np.concatenate(([0], np.log(dm.vp.z[1:]), np.log(dm.vp.values), np.log(dm.vs.values), np.log(dm.rh.values)))

    # perturbation step for each dimension
    dmodel = np.concatenate((dlogz  * np.ones_like(dm.vp.z),
                             dlogvp * np.ones_like(dm.vp.z),
                             dlogvs * np.ones_like(dm.vp.z),
                             dlogrh * np.ones_like(dm.vp.z),
                             ))
    #
    values0 = herrmanncaller.disperse(
        dm.vp.z, dm.vp.values, dm.vs.values, dm.rh.values)

    curves = []  # overwrite input !!
    for w, t, m, F, Iwtm in groupbywtm(waves, types, modes, freqs, np.arange(len(waves))):
        curve = Curve(wave=w, type=t, mode=m, freqs=F, values=values0[Iwtm], )
        curves.append(curve)
    logvalues0 = lognofail(values0)

    IZ  = np.arange(0*nlayer, 1*nlayer)
    IVP = np.arange(1*nlayer, 2*nlayer)
    IVS = np.arange(2*nlayer, 3*nlayer)
    IRH = np.arange(3*nlayer, 4*nlayer)

    assert (dmodel[IZ]  == dlogz).all()
    assert (dmodel[IVP] == dlogvp).all()
    assert (dmodel[IVS] == dlogvs).all()
    assert (dmodel[IRH] == dlogrh).all()

    # allocate empty matrix for all kernels on top of each other,
    # all dispersion points in columns
    DVADP = np.zeros((4 * nlayer, len(waves)), float) * np.nan

    # ----
    # parallel
    # ----
    def fun(i, modeli):
        logztopi, logvpi, logvsi, logrhi = \
            modeli[IZ], modeli[IVP], modeli[IVS], modeli[IRH]
        n = len(logztopi)
        ztopi = np.exp(logztopi)
        ztopi[0] = 0.

        ilayer = i % n
        H = dm.vp.z[-1] - dm.vp.z[0]
        if ilayer == n-1:
            Hi = 1.e50  # thickness of the half-space
        else:
            Hi = ztopi[ilayer + 1] - ztopi[ilayer]

        try:
            logvaluesi = lognofail(herrmanncaller.disperse(
                ztopi,
                np.exp(logvpi),
                np.exp(logvsi),
                np.exp(logrhi)))
        except CPiSDomainError as err:
            print("error during gradient computation %s" % str(err))
            return i, None

        # absolute sensitivity regardless the thickness differences
        DVAVPi = (logvaluesi - logvalues0) / (modeli[i] - model0[i])
        if norm:
            # sensitivity corrected from the layer thicknesses
            DVAVPi *= H / Hi

        return i, DVAVPi

    def gen():
        # loop over layers, then parameters. skip depth of surface
        for i in range(1, 4 * nlayer):
            modeli = model0.copy()
            modeli[i] += dmodel[i]
            yield Job(i, modeli)

    with MapSync(fun, gen()) as ma:
        for _, (i, DVAVPi), _, _ in ma:
            if DVAVPi is None:
                continue
            DVADP[i, :] = DVAVPi



    for w, t, m, F, Iwtm in groupbywtm(waves, types, modes, freqs, np.arange(len(waves))):
        # curve = Curve(wave=w, type=t, mode=m, freqs=F, values=np.exp(logvalues0[Iwtm]), )
        for curve in curves:
            if curve.wave == w and curve.type == t and curve.mode == m:
                break
        else:
            raise

        # recover the part of the full kernel corresponding to this specific curve and parameter
        dvadz  = DVADP[IZ,  :][:, Iwtm]
        dvadvp = DVADP[IVP, :][:, Iwtm]
        dvadvs = DVADP[IVS, :][:, Iwtm]
        dvadrh = DVADP[IRH, :][:, Iwtm]

        dvadz, dvadvp, dvadvs, dvadrh = \
            [np.ma.masked_where(np.isnan(_), _) for _ in
             [dvadz, dvadvp, dvadvs, dvadrh]]

        # # dlnc = dronlnc / dronlnvp dlnvp + drondlnc / dronlnvs dlnvs
        # # lnvp = lnvs + ln(vp/vs) => dlnvp / dln(vp/vs) = 1
        # # lnvs = lnvp - ln(vp/vs) => dlnvs / dln(vp/vs) = -1
        # # => dlnc / dln(vp/vs) = (dronlnc / dronlnvp) - (drondlnc / dronlnvs)
        # DLOGVADLOGPR = (DLOGVADLOGVP - DLOGVADLOGVS)

        if not relative:
            # dc/dvs = c/vs * dlnc/dlnvs
            dvadz[1:, :] *= curve.values / dm.vp.z[1:, None]
            dvadvp *= curve.values / dm.vp.values[:, None]
            dvadvs *= curve.values / dm.vs.values[:, None]
            dvadrh *= curve.values / dm.rh.values[:, None]

        # yield w, t, m, F, DLOGVADLOGZ, DLOGVADLOGVP, DLOGVADLOGVS, DLOGVADLOGRH # , DLOGVADLOGPR
        kernel_ztop = SensitivityKernel(parameter_name="ztop", curve=curve, dm=dm, kernel=dvadz,  all_curves=curves, norm=norm, relative=relative)
        kernel_vp   = SensitivityKernel(parameter_name="vp",   curve=curve, dm=dm, kernel=dvadvp, all_curves=curves, norm=norm, relative=relative)
        kernel_vs   = SensitivityKernel(parameter_name="vs",   curve=curve, dm=dm, kernel=dvadvs, all_curves=curves, norm=norm, relative=relative)
        kernel_rh   = SensitivityKernel(parameter_name="rh",   curve=curve, dm=dm, kernel=dvadrh, all_curves=curves, norm=norm, relative=relative)

        yield w, t, m, (kernel_ztop, kernel_vp, kernel_vs, kernel_rh, )

#
# def sker17_1(ztop, vp, vs, rh,
#     Waves, Types, Modes, Freqs, **kwargs):
#     """sker17_1 : same as sker17 with slightely more convenient input (no need to repeat wave, type and mode)
#
#     Waves is like ['L', 'L', 'R']
#     Types is like ['C', 'C', 'U']
#     Modes is like [ 0,   1,   0 ]
#     Freqs is like [fLC0, fLC1, fRU0] where f??? are frequency numpy arrays or lists
#
#     see sker17 for detailed input and output arguments
#     """
#     waves, types, modes, freqs = igroupbywtm(Waves, Types, Modes, Freqs)
#     for tup in sker17(ztop, vp, vs, rh,
#             waves, types, modes, freqs, **kwargs):
#         yield tup
#
#
# def sker17_2(depthmodel: depthmodel, curves: List[Curve], **kwargs):
#     """sker17_2 : same as sker17 with same arguments as HerrmannCaller
#
#     see sker17 for detailed input and output arguments
#     """
#
#     ztop = depthmodel.vp.ztop()
#     vp = depthmodel.vp.values
#     vs = depthmodel.vs.values
#     rh = depthmodel.rh.values
#
#     Waves, Types, Modes, Freqs = zip(*[(curve.wave, curve.type, curve.mode, curve.freqs) for curve in curves])
#     waves, types, modes, freqs = igroupbywtm(Waves, Types, Modes, Freqs)
#
#     for tup in sker17(ztop, vp, vs, rh,
#             waves, types, modes, freqs,
#             **kwargs):
#         yield tup


def main():

    help = '''sker17
    -m96          depthmodel to read 
    -RU0          rayleigh, group, mode 0 : expects 4 frequency arguments : fstart, fend, nfreq, fscale
    -RU1          rayleigh, group, mode 1 : expects 4 frequency arguments : fstart, fend, nfreq, fscale
    -RC0          rayleigh, phase, mode 0 : expects 4 frequency arguments : fstart, fend, nfreq, fscale          
    -LC0          love,     phase, mode 0 : expects 4 frequency arguments : fstart, fend, nfreq, fscale
    ...
    -norm         if mentionned the kernels are divided by the layer thickness
                  use it for depth models with irregular thicknesses 
    -png          save figures as pngfiles (overwrite if exists)
                  sker17_depthdisp.png
                  sker17_RU0_fstart_fend_nfreq_fscale.png
                  ...
    '''
    raise Exception('TODO : upgrade')
    if len(sys.argv) == 1:
        print(help)
        sys.exit()

    argv = readargv()
    # -----------------------------------:
    dm = depthmodel_from_mod96(argv['m96'][0])
    ztop = dm.vp.ztop()
    vp = dm.vp.values
    vs = dm.vs.values
    rh = dm.rh.values

    norm = "norm" in argv.keys()
    png = "png" in argv.keys()
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

    # ==== compute the dispersion curves
    hc = HerrmannCallerFromGroupedLists(Waves, Types, Modes, Freqs, h=0.005, ddc=0.005)
    with Timer('dispersion'):
        curves_out = hc(ztop, vp, vs, rh)

    # ==== display
    fig1 = plt.figure(figsize=(8,8))
    fig1.subplots_adjust(wspace=0.3)

    # depth model
    ax1 = fig1.add_subplot(223)
    dm.show(ax1)
    ax1.grid(True, linestyle=":", color="k")
    plt.legend()

    # disp curve
    ax2 = fig1.add_subplot(222)
    for curve in curves_out:
        # ax2.loglog(1. / fs, us, '+-', label="%s%s%d" % (w, t, m))
        curve.plot(ax2, "+-")
    ax2.set_ylabel('velocity (km/s)')
    ax2.grid(True, which="major")
    ax2.grid(True, which="minor")
    logtick(ax2, "xy")
    plt.legend()

    # ## sensitivity kernels
    if not png:
        fig1.show()

    sker_gen = sker17_1(ztop, vp, vs, rh,
                     Waves, Types, Modes, Freqs,
                     dz=0.001, dlogvs=.01,
                     dlogpr=.01, dlogrh=.01, norm=norm,
                     h=0.005, ddc=0.005)
                     
    for w, t, m, F, DLOGVADZ, DLOGVADLOGVS, DLOGVADLOGPR, DLOGVADLOGRH in sker_gen:
            
        # ------
        z_edges = np.hstack((ztop, [1.1 * ztop[-1]]))
        z_mid = np.hstack((0.5 * (ztop[1:] + ztop[:-1]), [1.05 * ztop[-1]]))
        F_edges = np.hstack((F[0] * 0.95, np.sqrt(F[1:] * F[:-1]), F[-1] * 1.05))

        # ------
        #vmax = abs(DLOGVADLOGVS).max()
        # #np.max([abs(DLOGVADZ).max(), abs(DLOGVADLOGVS).max(),
        # abs(DLOGVADLOGPR).max(), abs(DLOGVADLOGRH).max()])

        if not norm:
            # mask half space because it integrates the sensitivity over very thick layer => overestimated sensitivity
            for _ in DLOGVADZ, DLOGVADLOGVS, DLOGVADLOGPR, DLOGVADLOGRH:
                _[-1, :] = np.nan
            DLOGVADZ, DLOGVADLOGVS, DLOGVADLOGPR, DLOGVADLOGRH = \
                [np.ma.masked_where(np.isnan(_), _) \
                 for _ in [DLOGVADZ, DLOGVADLOGVS, DLOGVADLOGPR, DLOGVADLOGRH]]

        cmap = tomocmap1(w=0.01, W=0.2)  # cccfcmap3() #plt.cm.RdBu
        for M, p, q in zip([DLOGVADZ, DLOGVADLOGVS, DLOGVADLOGPR, DLOGVADLOGRH],
                           ["Z^{top}_i", "ln Vs_i", "ln (Vp/Vs)_i", r"ln \rho _i"],
                           ["Ztop", "lnVs", "lnVpaVs", "lnrho"]):
            ax3 = fig1.add_subplot(224, sharex=ax2, sharey=ax1)
            ax3.set_title('%s%s%d' % (w, t, m))

            vmax = abs(M).max()
            vmin = -vmax

            coll = plt.pcolormesh(1. / F_edges, z_edges, M,
                                  vmin=vmin, vmax=vmax, cmap=cmap)

            if M.max() - M.min():
                levels = np.logspace(-1., 2, 10)
                levels = np.hstack((-levels[::-1], [0], levels))
                plt.contour(1. / F, z_mid, M,
                            levels=levels,
                            colors="k")

            cax = fig1.add_axes((.91, .2, .01, .2))
            plt.colorbar(coll, cax=cax)

            ax3.set_xlabel('period (s)')
            ax3.set_xlim(minmax(1. / F_edges))
            ax3.set_ylim(minmax(z_edges))
            ax3.set_xscale('log')

            if norm:
                textonly(ax3, txt=r'$ \frac{H}{H_i} \, \frac{d ln%s_j}{d %s} $' % (t, p), loc=3, fontsize=16)
            else:
                textonly(ax3, txt=r'$ \frac{d ln%s_j}{d %s} $' % (t, p), loc=3, fontsize=16)

            if not ax3.yaxis_inverted():
                ax3.invert_yaxis()
            logtick(ax3, "x")

            # ------
            # plt.setp(ax1.get_xticklabels(), visible=False)
            # plt.setp(ax2.get_xticklabels(), visible=False)
            # plt.setp(ax2.get_yticklabels(), visible=False)
            # plt.setp(ax4.get_yticklabels(), visible=False)
            fig1.canvas.draw()

            if "png" in argv.keys():
                k = "%s%s%d" % (w, t, m)
                fout = 'sker17_%s_%s_%s_%s_%s_%s%s.png' % (k, argv[k][0], argv[k][1], argv[k][2], argv[k][3], q, "_norm" if norm else "")
                print(fout)
                fig1.savefig(fout, dpi=300)
            else:
                input('pause : press enter to plot the next wave type and mode')
            cax.cla()
            ax3.cla()

    # --------------------
    if "png" not in argv.keys():
        input('bye')


# -----------------------------
if __name__ == "__main__":
    main()
