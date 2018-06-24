from srfpython.standalone.display import plt, gcf, gca, pause, showme, Ntick, logtick, makecolorbar
from srfpython.depthdisp.dispcurves import surf96reader, mklaws
from srfpython.depthdisp.depthmodels import depthmodel1D, depthmodel_from_mod96
import numpy as np


# _____________________________________________
class DepthDispDisplay(object):
    """convenient tool to display depthmodels and dispersion curves on the same figure"""

    def __init__(self, fig=None, targetfile=None):
        if fig is None:
            self.fig = plt.figure(figsize=(12, 6))#figsize=(18, 10))
            self.fig.subplots_adjust(wspace=0.05)
        else:
            self.fig = fig

        self.axdepth = {}
        self.axdepth['VP'] = self.fig.add_subplot(1, 5, 1, title="$V_P\,(km/s)$", ylabel="depth (km)")
        self.axdepth['VS'] = self.fig.add_subplot(1, 5, 2, title="$V_S\,(km/s)$",sharey=self.axdepth['VP'])
        self.axdepth['PR'] = self.fig.add_subplot(1, 5, 3, title=r"$V_P/V_S$", sharey=self.axdepth['VP'])
        self.axdepth['RH'] = self.fig.add_subplot(1, 5, 4, title=r"$\rho\,(g/cm^3)$", sharey=self.axdepth['VP'])

        if targetfile is None:
            self.axdisp = {}
            self.axdisp['RU0'] = self.fig.add_subplot(4, 5, 5, title="RU0", ylabel="grpvel (km/s)")
            self.axdisp['RU1'] = self.fig.add_subplot(4, 5, 10, title="RU1", ylabel="grpvel (km/s)", sharex=self.axdisp['RU0'], sharey=self.axdisp['RU0'])
            self.axdisp['RC0'] = self.fig.add_subplot(4, 5, 15, title="RC0", ylabel="phsvel (km/s)", sharex=self.axdisp['RU0'], sharey=self.axdisp['RU0'])
            self.axdisp['RC1'] = self.fig.add_subplot(4, 5, 20, title="RC1", ylabel="phsvel (km/s)", sharex=self.axdisp['RU0'], sharey=self.axdisp['RU0'])
            for key, ax in self.axdisp.items():
                ax.loglog([1.],[1.]) #fuck matplotlib v2
                ax.yaxis.set_label_position("right")
                ax.yaxis.tick_right()
        else:
            #adjust the suplots according to the target data
            s = surf96reader(targetfile)
            Ndisp = max([4, len(list(s.wtm()))])

            self.axdisp = {}
            share = None
            for n, (w, t, m)  in enumerate(s.wtm()):
                ax = self.fig.add_subplot(Ndisp, 5, (n+1)*5,
                                          title="%s%s%d" % (w, t, m),
                                          sharex=share, sharey=share,
                                          #ylabel="velocity (km/s)")
                                          ylabel="%s (km/s)" % (["grpvel", "phsvel"][int(t=="C")]))
                ax.yaxis.set_label_position("right")
                ax.yaxis.tick_right()

                self.axdisp["%s%s%d" % (w.upper(), t.upper(), m)] = share = ax
                plt.setp(ax.get_xticklabels(), visible = False)
            plt.setp(ax.get_xticklabels(), visible=True)
            ax.set_xlabel('$period\,(s)$')

        pos = self.axdisp[self.axdisp.keys()[-1]].get_position()
        # self.cax = self.fig.add_axes((pos.x0, 0.12, pos.width, 0.01))
        self.cax = self.fig.add_axes((pos.x0, pos.y0 - 0.12, pos.width, 0.01))

        plt.setp(self.axdepth['VS'].get_yticklabels(), visible=False)
        plt.setp(self.axdepth['PR'].get_yticklabels(), visible=False)
        plt.setp(self.axdepth['RH'].get_yticklabels(), visible=False)
        self.axconv = None #Not needed here

        if not self.axdepth['VP'].yaxis_inverted():
            self.axdepth['VP'].invert_yaxis()

        # initiate collections data
        # self.clear_collections()

    def clear_collections(self):
        for k in self.axdepth.keys():
            self.deptcoll[k] = []
        for k in self.axdisp.keys():
            self.dispcoll[k] = []

    def colorbar(self, vmin, vmax, cmap, **kwargs):
        cb = makecolorbar(vmin=vmin, vmax=vmax, cmap=cmap)
        self.fig.colorbar(cb, cax=self.cax, **kwargs)

    def cla(self):
        for _, ax in self.axdisp.items(): ax.cla()
        for _, ax in self.axdepth.items(): ax.cla()
        self.axconv.cla()

    def cliplim(self):
        vplim = self.axdepth['VP'].get_xlim()
        vslim = self.axdepth['VS'].get_xlim()
        prlim = self.axdepth['PR'].get_xlim()
        rhlim = self.axdepth['RH'].get_xlim()
        vplim = np.clip(vplim, 0., np.inf)
        vslim = np.clip(vslim, 0., np.inf)
        prlim = np.clip(prlim, 1., np.inf)
        rhlim = np.clip(rhlim, 1., np.inf)
        self.axdepth['VP'].set_xlim(vplim)
        self.axdepth['VS'].set_xlim(vslim)
        self.axdepth['PR'].set_xlim(prlim)
        self.axdepth['RH'].set_xlim(rhlim)

    def grid(self):
        for _, ax in self.axdisp.items():  ax.grid(True, linestyle=":")
        for _, ax in self.axdepth.items(): ax.grid(True, linestyle=":")
        if self.axconv is not None:
            self.axconv.grid(True, linestyle=":")

    def tick(self):
        for _, ax in self.axdisp.items():  logtick(ax, "xy")
        for _, ax in self.axdepth.items(): Ntick(ax, 4, "x")

    def plotmodel(self, ztop, vp, vs, rh, color="k", alpha=0.2, showvp=True, showvs=True, showrh=True, showpr=True,
                  **kwargs):
        if showpr: pr = depthmodel1D(ztop, vp / vs)
        if showvp: vp = depthmodel1D(ztop, vp)
        if showvs: vs = depthmodel1D(ztop, vs)
        if showrh: rh = depthmodel1D(ztop, rh)

        if showpr: pr.show(self.axdepth['PR'], color=color, alpha=alpha, **kwargs)
        if showvp: vp.show(self.axdepth['VP'], color=color, alpha=alpha, **kwargs)
        if showvs: vs.show(self.axdepth['VS'], color=color, alpha=alpha, **kwargs)
        if showrh: rh.show(self.axdepth['RH'], color=color, alpha=alpha, **kwargs)

    def plotm96(self, m96, **kwargs):
        dm = depthmodel_from_mod96(m96)
        ztop, vp, vs, rh = dm.ztopvpvsrh()
        self.plotmodel(ztop, vp, vs, rh, **kwargs)

    def plotdisp(self, waves, types, modes, freqs, values, dvalues=None, color="k", alpha=0.2, **kwargs):
        for law in mklaws(waves, types, modes, freqs, values, dvalues=dvalues):
            ax = self.axdisp["%s%s%d" % (law.wave.upper(), law.type.upper(), law.mode)]
            law.show(ax, period=True, showdvalue=True, color=color, alpha=alpha, **kwargs)

    def plotvspdf(self, vspdf, **kwargs):
        vspdf.show(self.axdepth['VS'], **kwargs)

    def plotvppdf(self, vppdf, **kwargs):
        vppdf.show(self.axdepth['VP'], **kwargs)

    def plotrhpdf(self, rhpdf, **kwargs):
        rhpdf.show(self.axdepth['RH'], **kwargs)

    def plotprpdf(self, prpdf, **kwargs):
        prpdf.show(self.axdepth['PR'], **kwargs)

    def plotru0pdf(self, ru0pdf, **kwargs):
        ru0pdf.show(self.axdisp['RU0'], **kwargs)

    def plotru1pdf(self, ru1pdf, **kwargs):
        ru1pdf.show(self.axru1, **kwargs)

    def plots96(self, s96, **kwargs):
        s = surf96reader(s96)
        waves, types, modes, freqs, values, dvalues = s.wtmfvd()
        self.plotdisp(waves, types, modes, freqs, values, dvalues=dvalues, **kwargs)

    def set_plim(self, plim):
        for _, ax in self.axdisp.items():
            ax.set_xlim(plim)

    def set_vlim(self, vlim):
        for _, ax in self.axdisp.items():
            ax.set_ylim(vlim)

    def set_vplim(self, lim):
        self.axdepth['VP'].set_xlim(lim)

    def set_vslim(self, lim):
        self.axdepth['VS'].set_xlim(lim)

    def set_prlim(self, lim):
        self.axdepth['PR'].set_xlim(lim)

    def set_rhlim(self, lim):
        self.axdepth['RH'].set_xlim(lim)

    def set_zlim(self, zlim):
        for _, ax in self.axdepth.items():
            ax.set_ylim(max(abs(zlim)), min(abs(zlim)))

# _____________________________________________
class DepthDispDisplayCompact(DepthDispDisplay):
    """display only vs and the dispersion curves"""
    def __init__(self, fig=None, targetfile=None):
        if fig is None:
            self.fig = plt.figure(figsize=(5, 6))#figsize=(18, 10))
            self.fig.subplots_adjust(wspace=0.05)
        else:
            self.fig = fig

        self.axdepth['VS'] = self.fig.add_subplot(1, 2, 1, title="$V_S\,(km/s)$", ylabel="depth (km)")
        self.axdepth['VP'] = self.axdepth['PR'] = self.axdepth['RH'] = None

        if targetfile is None:
            self.axdisp = {}
            self.axdisp['RU0'] = self.fig.add_subplot(3, 2, 2, title="RU0", ylabel="grpvel (km/s)")
            self.axdisp['RU1'] = self.fig.add_subplot(3, 2, 4, title="RU1", ylabel="grpvel (km/s)", sharex=self.axdisp['RU0'], sharey=self.axdisp['RU0'])

            for _, ax in self.axdisp.items():
                ax.loglog([1.],[1.]) #fuck matplotlib v2
                ax.yaxis.set_label_position("right")
                ax.yaxis.tick_right()
        else:
            #adjust the suplots according to the target data
            s = surf96reader(targetfile)
            Ndisp = max([3, len(list(s.wtm()))])

            self.axdisp = {}
            share = None
            for n, (w, t, m)  in enumerate(s.wtm()):
                ax = self.fig.add_subplot(Ndisp, 2, (n+1)*2,
                                          title="%s%s%d" % (w, t, m),
                                          sharex=share, sharey=share,
                                          #ylabel="velocity (km/s)")
                                          ylabel="%s (km/s)" % (["grpvel", "phsvel"][int(t=="C")]))
                ax.yaxis.set_label_position("right")
                ax.yaxis.tick_right()
                share=ax
                self.axdisp["%s%s%d" % (w.upper(), t.upper(), m)] = share = ax
                plt.setp(ax.get_xticklabels(), visible = False)
            plt.setp(ax.get_xticklabels(), visible=True)
            ax.set_xlabel('$period\,(s)$')

        pos = self.axdisp[self.axdisp.keys()[-1]].get_position()
        #self.cax = self.fig.add_axes((pos.x0, 0.12, pos.width, 0.01))
        self.cax = self.fig.add_axes((pos.x0, pos.y0 - 0.12, pos.width, 0.01))

        # plt.setp(self.axdepth['VS'].get_yticklabels(), visible=False)
        self.axconv = None #Not needed here

        if not self.axdepth['VS'].yaxis_inverted():
            self.axdepth['VS'].invert_yaxis()

    # _____________________________________________
    def plotmodel(self, *args, **kwargs):
        kwargs.setdefault('showpr', False)
        kwargs.setdefault('showvp', False)
        kwargs.setdefault('showrh', False)

        DepthDispDisplay.plotmodel(self, *args, **kwargs)
