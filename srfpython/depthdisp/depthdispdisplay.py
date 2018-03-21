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

        self.axvp = self.fig.add_subplot(1, 5, 1, title="$V_P\,(km/s)$", ylabel="depth (km)")
        self.axvs = self.fig.add_subplot(1, 5, 2, title="$V_S\,(km/s)$", sharey=self.axvp)  # , sharex = self.axvp)
        self.axpr = self.fig.add_subplot(1, 5, 3, title=r"$V_P/V_S$", sharey=self.axvp)
        self.axrh = self.fig.add_subplot(1, 5, 4, title=r"$\rho\,(g/cm^3)$", sharey=self.axvp)

        if targetfile is None:
            self.axru0 = self.fig.add_subplot(4, 5, 5,  title="RU0", ylabel="grpvel (km/s)")
            self.axru1 = self.fig.add_subplot(4, 5, 10, title="RU1", ylabel="grpvel (km/s)", sharex=self.axru0, sharey=self.axru0)
            self.axrc0 = self.fig.add_subplot(4, 5, 15, title="RC0", ylabel="phsvel (km/s)", sharex=self.axru0, sharey=self.axru0)
            self.axrc1 = self.fig.add_subplot(4, 5, 20, title="RC1", ylabel="phsvel (km/s)", sharex=self.axru0, sharey=self.axru0)
            self.axdisp = [self.axru0, self.axru1, self.axrc0, self.axrc1]
            for ax in self.axdisp:
                ax.loglog([1.],[1.]) #fuck matplotlib v2
                ax.yaxis.set_label_position("right")
                ax.yaxis.tick_right()
        else:
            #adjust the suplots according to the target data
            s = surf96reader(targetfile)
            Ndisp = max([4, len(list(s.wtm()))])

            self.axdisp = []
            share = None
            for n, (w, t, m)  in enumerate(s.wtm()):
                ax = self.fig.add_subplot(Ndisp, 5, (n+1)*5,
                                          title="%s%s%d" % (w, t, m),
                                          sharex=share, sharey=share,
                                          #ylabel="velocity (km/s)")
                                          ylabel="%s (km/s)" % (["grpvel", "phsvel"][int(t=="C")]))
                ax.yaxis.set_label_position("right")
                ax.yaxis.tick_right()
                share=ax
                self.__setattr__("ax%s%s%d" % (w.lower(), t.lower(), m), ax)
                self.axdisp.append(eval("self.ax%s%s%d" % (w.lower(), t.lower(), m)))
                plt.setp(ax.get_xticklabels(), visible = False)
            plt.setp(ax.get_xticklabels(), visible=True)
            ax.set_xlabel('$period\,(s)$')

        pos = self.axdisp[-1].get_position()
        #self.cax = self.fig.add_axes((pos.x0, 0.12, pos.width, 0.01))
        self.cax = self.fig.add_axes((pos.x0, pos.y0 - 0.12, pos.width, 0.01))

        plt.setp(self.axvs.get_yticklabels(), visible=False)
        plt.setp(self.axpr.get_yticklabels(), visible=False)
        plt.setp(self.axrh.get_yticklabels(), visible=False)
        self.axconv = None #Not needed here

        self.axdepth = [self.axvp, self.axvs, self.axpr, self.axrh]
        if not self.axvp.yaxis_inverted():
            self.axvp.invert_yaxis()

    def colorbar(self, vmin, vmax, cmap, **kwargs):
        cb = makecolorbar(vmin=vmin, vmax=vmax, cmap=cmap)
        self.fig.colorbar(cb, cax=self.cax, **kwargs)

    def cla(self):
        for ax in self.axdisp: ax.cla()
        for ax in self.axdepth: ax.cla()
        self.axconv.cla()

    def cliplim(self):
        vplim = self.axvp.get_xlim()
        vslim = self.axvs.get_xlim()
        prlim = self.axpr.get_xlim()
        rhlim = self.axrh.get_xlim()
        vplim = np.clip(vplim, 0., np.inf)
        vslim = np.clip(vslim, 0., np.inf)
        prlim = np.clip(prlim, 1., np.inf)
        rhlim = np.clip(rhlim, 1., np.inf)
        self.axvp.set_xlim(vplim)
        self.axvs.set_xlim(vslim)
        self.axpr.set_xlim(prlim)
        self.axrh.set_xlim(rhlim)

    def grid(self):
        for ax in self.axdisp:  ax.grid(True, linestyle=":")
        for ax in self.axdepth: ax.grid(True, linestyle=":")
        if self.axconv is not None:
            self.axconv.grid(True, linestyle=":")

    def tick(self):
        for ax in self.axdisp:  logtick(ax, "xy")
        for ax in self.axdepth: Ntick(ax, 4, "x")

    def plotmodel(self, ztop, vp, vs, rh, color="k", alpha=0.2, showvp=True, showvs=True, showrh=True, showpr=True,
                  **kwargs):
        if showpr: pr = depthmodel1D(ztop, vp / vs)
        if showvp: vp = depthmodel1D(ztop, vp)
        if showvs: vs = depthmodel1D(ztop, vs)
        if showrh: rh = depthmodel1D(ztop, rh)

        if showpr: pr.show(self.axpr, color=color, alpha=alpha, **kwargs)
        if showvp: vp.show(self.axvp, color=color, alpha=alpha, **kwargs)
        if showvs: vs.show(self.axvs, color=color, alpha=alpha, **kwargs)
        if showrh: rh.show(self.axrh, color=color, alpha=alpha, **kwargs)

    def plotm96(self, m96, **kwargs):
        dm = depthmodel_from_mod96(m96)
        ztop, vp, vs, rh = dm.ztopvpvsrh()
        self.plotmodel(ztop, vp, vs, rh, **kwargs)

    def plotdisp(self, waves, types, modes, freqs, values, dvalues=None, color="k", alpha=0.2, **kwargs):
        for law in mklaws(waves, types, modes, freqs, values, dvalues=dvalues):
            ax = eval("self.ax%s%s%d" % (law.wave.lower(), law.type.lower(), law.mode))
            law.show(ax, period=True, showdvalue=True, color=color, alpha=alpha, **kwargs)

    def plotvspdf(self, vspdf, **kwargs):
        vspdf.show(self.axvs, **kwargs)

    def plotvppdf(self, vppdf, **kwargs):
        vppdf.show(self.axvp, **kwargs)

    def plotrhpdf(self, rhpdf, **kwargs):
        rhpdf.show(self.axrh, **kwargs)

    def plotprpdf(self, prpdf, **kwargs):
        prpdf.show(self.axpr, **kwargs)

    def plotru0pdf(self, ru0pdf, **kwargs):
        ru0pdf.show(self.axru0, **kwargs)

    def plotru1pdf(self, ru1pdf, **kwargs):
        ru1pdf.show(self.axru1, **kwargs)

    def plots96(self, s96, **kwargs):
        s = surf96reader(s96)
        waves, types, modes, freqs, values, dvalues = s.wtmfvd()
        self.plotdisp(waves, types, modes, freqs, values, dvalues=dvalues, **kwargs)

    def set_plim(self, plim):
        for ax in self.axdisp:
            ax.set_xlim(plim)

    def set_vlim(self, vlim):
        for ax in self.axdisp:
            ax.set_ylim(vlim)

    def set_vplim(self, lim):
        self.axvp.set_xlim(lim)

    def set_vslim(self, lim):
        self.axvs.set_xlim(lim)

    def set_prlim(self, lim):
        self.axpr.set_xlim(lim)

    def set_rhlim(self, lim):
        self.axrh.set_xlim(lim)

    def set_zlim(self, zlim):
        for ax in self.axdepth:
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

        self.axvs = self.fig.add_subplot(1, 2, 1, title="$V_S\,(km/s)$", ylabel="depth (km)")
        self.axvp = self.axpr = self.axrh = None

        if targetfile is None:
            self.axru0 = self.fig.add_subplot(3, 2, 2, title="RU0", ylabel="grpvel (km/s)")
            self.axru1 = self.fig.add_subplot(3, 2, 4, title="RU1", ylabel="grpvel (km/s)", sharex=self.axru0, sharey=self.axru0)
            # self.axrc0 = self.fig.add_subplot(4, 2, 6, title="RC0", ylabel="phsvel (km/s)", sharex=self.axru0, sharey=self.axru0)
            # self.axrc1 = self.fig.add_subplot(4, 2, 8, title="RC1", ylabel="phsvel (km/s)", sharex=self.axru0, sharey=self.axru0)
            self.axdisp = [self.axru0, self.axru1] #, self.axrc0, self.axrc1]
            for ax in self.axdisp:
                ax.loglog([1.],[1.]) #fuck matplotlib v2
                ax.yaxis.set_label_position("right")
                ax.yaxis.tick_right()
        else:
            #adjust the suplots according to the target data
            s = surf96reader(targetfile)
            Ndisp = max([3, len(list(s.wtm()))])

            self.axdisp = []
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
                self.__setattr__("ax%s%s%d" % (w.lower(), t.lower(), m), ax)
                self.axdisp.append(eval("self.ax%s%s%d" % (w.lower(), t.lower(), m)))
                plt.setp(ax.get_xticklabels(), visible = False)
            plt.setp(ax.get_xticklabels(), visible=True)
            ax.set_xlabel('$period\,(s)$')

        pos = self.axdisp[-1].get_position()
        #self.cax = self.fig.add_axes((pos.x0, 0.12, pos.width, 0.01))
        self.cax = self.fig.add_axes((pos.x0, pos.y0 - 0.12, pos.width, 0.01))

        # plt.setp(self.axvs.get_yticklabels(), visible=False)
        self.axconv = None #Not needed here

        self.axdepth = [self.axvs]
        if not self.axvs.yaxis_inverted():
            self.axvs.invert_yaxis()

    # _____________________________________________
    def plotmodel(self, *args, **kwargs):
        kwargs.setdefault('showpr', False)
        kwargs.setdefault('showvp', False)
        kwargs.setdefault('showrh', False)

        DepthDispDisplay.plotmodel(self, *args, **kwargs)
