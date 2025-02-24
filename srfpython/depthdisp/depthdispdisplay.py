from srfpython.standalone.display import plt, gcf, gca, pause, showme, Ntick, logtick, makecolorbar
from srfpython.depthdisp.dispcurves import surf96reader, mklaws
from srfpython.depthdisp.depthmodels import depthmodel1D, depthmodel_from_mod96
from matplotlib.collections import LineCollection
import numpy as np


class DepthDispDisplay(object):
    """Tool to display depthmodels and dispersion curves on the same figure"""
    dist_unit = "km"
    time_unit = "s"
    density_unit = "g/cm3"
    dist_scales = {"km": 1.0, "m": 1000.}
    time_scales = {"s": 1.0, "ms": 1000.}
    density_scales = {"g/cm3": 1.0, "kg/m3": 1000.}
    period = True

    def __init__(self, fig=None, targetfile=None):
        if fig is None:
            self.fig = plt.figure(figsize=(12, 6))#figsize=(18, 10))
            self.fig.subplots_adjust(wspace=0.05)
        else:
            fig.clf()
            self.fig = fig

        self.axdepth = {}
        self.axdepth['VP'] = self.fig.add_subplot(1, 5, 1, title="$V_P\,(%s/%s)$" % (self.dist_unit, self.time_unit), ylabel="depth (%s)" % self.dist_unit)
        self.axdepth['VS'] = self.fig.add_subplot(1, 5, 2, title="$V_S\,(%s/%s)$" % (self.dist_unit, self.time_unit), sharey=self.axdepth['VP'])
        self.axdepth['PR'] = self.fig.add_subplot(1, 5, 3, title=r"$V_P/V_S$", sharey=self.axdepth['VP'])
        self.axdepth['RH'] = self.fig.add_subplot(1, 5, 4, title=r"$\rho\,(%s)$" % self.density_unit, sharey=self.axdepth['VP'])

        if targetfile is None:
            self.axdisp = {}
            self.axdisp['RU0'] = self.fig.add_subplot(4, 5, 5, title="RU0", ylabel="grpvel (%s/%s)" % (self.dist_unit, self.time_unit))
            self.axdisp['RU1'] = self.fig.add_subplot(4, 5, 10, title="RU1", ylabel="grpvel (%s/%s)" % (self.dist_unit, self.time_unit), sharex=self.axdisp['RU0'], sharey=self.axdisp['RU0'])
            self.axdisp['RC0'] = self.fig.add_subplot(4, 5, 15, title="RC0", ylabel="phsvel (%s/%s)" % (self.dist_unit, self.time_unit), sharex=self.axdisp['RU0'], sharey=self.axdisp['RU0'])
            self.axdisp['RC1'] = self.fig.add_subplot(4, 5, 20, title="RC1", ylabel="phsvel (%s/%s)" % (self.dist_unit, self.time_unit), sharex=self.axdisp['RU0'], sharey=self.axdisp['RU0'])
            for key, ax in self.axdisp.items():
                ax.loglog([1.],[1.]) #fuck matplotlib v2
                ax.yaxis.set_label_position("right")
                ax.yaxis.tick_right()
        else:
            # adjust the suplots according to the target data
            s = surf96reader(targetfile)
            Ndisp = max([4, len(list(s.wtm()))])

            self.axdisp = {}
            share = None
            wtms = list(s.wtm())
            for n, (w, t, m) in enumerate(wtms):

                ax = self.fig.add_subplot(Ndisp, 5, (n+1)*5,
                                          title="%s%s%d" % (w, t, m),
                                          sharex=share, sharey=share,
                                          #ylabel="velocity (km/s)")
                                          ylabel="%s (%s/%s)" % (["grpvel", "phsvel"][int(t=="C")], self.dist_unit, self.time_unit))
                ax.yaxis.set_label_position("right")
                ax.yaxis.tick_right()

                self.axdisp["%s%s%d" % (w.upper(), t.upper(), m)] = share = ax
                if n < len(wtms)-1:
                    plt.setp(ax.get_xticklabels(), visible=False)

            if self.period:
                ax.set_xlabel('$Period\,(%s)$' % self.time_unit)
            else:
                ax.set_xlabel('$Frequency\,(%s^{-1})$' % self.time_unit)

        pos = ax.get_position()
        self.cax = self.fig.add_axes((pos.x0, np.max([0.12, pos.y0 - 0.12]), pos.width, 0.01))

        plt.setp(self.axdepth['VS'].get_yticklabels(), visible=False)
        plt.setp(self.axdepth['PR'].get_yticklabels(), visible=False)
        plt.setp(self.axdepth['RH'].get_yticklabels(), visible=False)
        self.axconv = None #Not needed here

        if not self.axdepth['VP'].yaxis_inverted():
            self.axdepth['VP'].invert_yaxis()

        # initiate collections data
        self.clear_collections()

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
        for _, ax in self.axdepth.items():
            if ax is not None:
                ax.grid(True, linestyle=":")
        if self.axconv is not None:
            self.axconv.grid(True, linestyle=":")

    def tick(self):
        for _, ax in self.axdisp.items():
            logtick(ax, "xy")

        for _, ax in self.axdepth.items():
            if ax is not None:
                Ntick(ax, 4, "x")

    def scale_and_plot_depthmodel1D(self, dm1d, which, **kwargs):

        dist_scale = self.dist_scales[self.dist_unit]
        time_scale = self.time_scales[self.time_unit]
        density_scale = self.density_scales[self.density_unit]

        scale = {"pr": 1.0,
                 "vp": dist_scale / time_scale,
                 "vs": dist_scale / time_scale,
                 "rh": density_scale}[which.lower()]

        # dm1d = dm1d.copy()
        dm1d.z *= dist_scale
        dm1d.values *= scale
        dm1d.show(
            ax=self.axdepth[which.upper()],
            **kwargs
            )

    def plot_depthmodel(self, dm, which, **kwargs):
        raise NotImplementedError

    def plotmodel(self, ztop, vp, vs, rh, color="k", alpha=0.2,
                  showvp=True, showvs=True, showrh=True, showpr=True,
                  **kwargs):

        if showpr:
            pr = depthmodel1D(ztop, vp / vs)
            self.scale_and_plot_depthmodel1D(pr, which="PR", color=color, alpha=alpha, **kwargs)

        if showvp:
            vp = depthmodel1D(ztop, vp)
            self.scale_and_plot_depthmodel1D(vp, which="VP", color=color, alpha=alpha, **kwargs)

        if showvs:
            vs = depthmodel1D(ztop, vs)
            self.scale_and_plot_depthmodel1D(vs, which="VS", color=color, alpha=alpha, **kwargs)

        if showrh:
            rh = depthmodel1D(ztop, rh)
            self.scale_and_plot_depthmodel1D(rh, which="RH", color=color, alpha=alpha, **kwargs)

    def plotm96(self, m96, **kwargs):
        dm = depthmodel_from_mod96(m96)

        ztop, vp, vs, rh = dm.ztopvpvsrh()
        # => plotmodel will apply scaling
        self.plotmodel(ztop, vp, vs, rh, **kwargs)

    def plotdisp(self, waves, types, modes, freqs, values, dvalues=None, color="k", alpha=0.2, **kwargs):

        dist_scale = self.dist_scales[self.dist_unit]
        time_scale = self.time_scales[self.time_unit]
        density_scale = self.density_scales[self.density_unit]

        for law in mklaws(waves, types, modes, freqs / time_scale, values * dist_scale / time_scale,
                          dvalues=dvalues if dvalues is None else dvalues * dist_scale / time_scale):
            ax = self.axdisp["%s%s%d" % (law.wave.upper(), law.type.upper(), law.mode)]
            law.show(ax, period=self.period, showdvalue=True, color=color, alpha=alpha, **kwargs)

    def plotvspdf(self, vspdf, **kwargs):
        # from srfpython.depthdisp.depthpdfs import depthpdf
        dist_scale = self.dist_scales[self.dist_unit]
        time_scale = self.time_scales[self.time_unit]
        density_scale = self.density_scales[self.density_unit]

        vspdf = vspdf.copy()
        vspdf.v *= dist_scale / time_scale
        vspdf.z *= dist_scale
        vspdf.show(self.axdepth['VS'], **kwargs)

    def plotvppdf(self, vppdf, **kwargs):
        dist_scale = self.dist_scales[self.dist_unit]
        time_scale = self.time_scales[self.time_unit]
        density_scale = self.density_scales[self.density_unit]

        vppdf = vppdf.copy()
        vppdf.v *= dist_scale / time_scale
        vppdf.z *= dist_scale

        vppdf.show(self.axdepth['VP'], **kwargs)

    def plotrhpdf(self, rhpdf, **kwargs):
        dist_scale = self.dist_scales[self.dist_unit]
        time_scale = self.time_scales[self.time_unit]
        density_scale = self.density_scales[self.density_unit]

        rhpdf = rhpdf.copy()
        rhpdf.v *= dist_scale / time_scale
        rhpdf.z *= dist_scale

        rhpdf.show(self.axdepth['RH'], **kwargs)

    def plotprpdf(self, prpdf, **kwargs):
        dist_scale = self.dist_scales[self.dist_unit]
        time_scale = self.time_scales[self.time_unit]
        density_scale = self.density_scales[self.density_unit]

        prpdf = prpdf.copy()
        prpdf.v *= dist_scale / time_scale
        prpdf.z *= dist_scale
        prpdf.show(self.axdepth['PR'], **kwargs)

    def plotru0pdf(self, ru0pdf, **kwargs):
        from srfpython.depthdisp.disppdfs import disppdf
        dist_scale = self.dist_scales[self.dist_unit]
        time_scale = self.time_scales[self.time_unit]
        density_scale = self.density_scales[self.density_unit]

        ru0pdf = ru0pdf.copy()
        ru0pdf.v *= dist_scale / time_scale
        ru0pdf.f *= 1. / time_scale
        ru0pdf.show(self.axdisp['RU0'], **kwargs)

    def plotru1pdf(self, ru1pdf, **kwargs):
        from srfpython.depthdisp.disppdfs import disppdf
        dist_scale = self.dist_scales[self.dist_unit]
        time_scale = self.time_scales[self.time_unit]
        density_scale = self.density_scales[self.density_unit]

        ru1pdf = ru1pdf.copy()
        ru1pdf.v *= dist_scale / time_scale
        ru1pdf.f *= 1. / time_scale
        ru1pdf.show(self.axru1, **kwargs)

    def plots96(self, s96, **kwargs):
        s = surf96reader(s96)

        waves, types, modes, freqs, values, dvalues = s.wtmfvd()
        # already scaled in plotdisp
        self.plotdisp(waves, types, modes, freqs, values, dvalues=dvalues, **kwargs)

    def set_plim(self, plim):
        time_scale = self.time_scales[self.time_unit]
        pmin, pmax = plim

        for _, ax in self.axdisp.items():
            if self.period:
                ax.set_xlim((pmin * time_scale, pmax * time_scale))
            else:
                ax.set_xlim((1. / time_scale / pmax, 1. / time_scale / pmin))

    def set_vlim(self, vlim):
        dist_scale = self.dist_scales[self.dist_unit]
        time_scale = self.time_scales[self.time_unit]
        vmin, vmax = vlim

        for _, ax in self.axdisp.items():
            ax.set_ylim((vmin * dist_scale / time_scale, vmax * dist_scale / time_scale))

    def set_vplim(self, lim):
        dist_scale = self.dist_scales[self.dist_unit]
        time_scale = self.time_scales[self.time_unit]
        vmin, vmax = lim

        self.axdepth['VP'].set_xlim(vmin * dist_scale / time_scale, vmax * dist_scale / time_scale)

    def set_vslim(self, lim):
        dist_scale = self.dist_scales[self.dist_unit]
        time_scale = self.time_scales[self.time_unit]
        vmin, vmax = lim
        self.axdepth['VS'].set_xlim(vmin * dist_scale / time_scale, vmax * dist_scale / time_scale)

    def set_prlim(self, lim):
        self.axdepth['PR'].set_xlim(lim)

    def set_rhlim(self, lim):
        densit_scale = self.dist_scales[self.density_unit]

        vmin, vmax = lim
        self.axdepth['RH'].set_xlim(vmin * densit_scale, vmax * densit_scale)

    def set_zlim(self, zlim):
        dist_scale = self.dist_scales[self.dist_unit]
        zmin = max(abs(zlim)) * dist_scale
        zmax = min(abs(zlim)) * dist_scale
        self.axdepth['VS'].set_ylim(zmin, zmax)

    def clear_collections(self):
        self.deptcoll = {}
        self.dispcoll = {}
        for k in self.axdepth.keys():
            self.deptcoll[k] = {'segments': [], 'colorvalues': []}
        for k in self.axdisp.keys():
            self.dispcoll[k] = {'segments': [], 'colorvalues': []}

    def addmodel(self, ztop, vp, vs, rh, colorvalue,
                 showvp=True, showvs=True, showrh=True, showpr=True):
        """same as plotmodel but the data are added to self.depthcoll for display as LineCollections
        :param self:
        :param ztop:
        :param vp:
        :param vs:
        :param rh:
        :param colorvalue:
        :param showvp:
        :param showvs:
        :param showrh:
        :param showpr:
        :return:
        """
        dist_scale = self.dist_scales[self.dist_unit]
        time_scale = self.time_scales[self.time_unit]
        density_scale = self.density_scales[self.density_unit]
        if showpr:
            pr = depthmodel1D(ztop * dist_scale, vp / vs)
            self.deptcoll['PR']['segments'].append(np.column_stack(pr.dontshow()))
            self.deptcoll['PR']['colorvalues'].append(colorvalue)
        if showvp:
            vp = depthmodel1D(ztop * dist_scale, vp * dist_scale / time_scale)
            self.deptcoll['VP']['segments'].append(np.column_stack(vp.dontshow()))
            self.deptcoll['VP']['colorvalues'].append(colorvalue)
        if showvs:
            vs = depthmodel1D(ztop * dist_scale, vs * dist_scale / time_scale)
            self.deptcoll['VS']['segments'].append(np.column_stack(vs.dontshow()))
            self.deptcoll['VS']['colorvalues'].append(colorvalue)
        if showrh:
            rh = depthmodel1D(ztop * dist_scale, rh * density_scale)
            self.deptcoll['RH']['segments'].append(np.column_stack(rh.dontshow()))
            self.deptcoll['RH']['colorvalues'].append(colorvalue)

    def adddisp(self, waves, types, modes, freqs, values, colorvalue):

        """add dispsersion curve to the right collection for massive display"""
        dist_scale = self.dist_scales[self.dist_unit]
        time_scale = self.time_scales[self.time_unit]
        density_scale = self.density_scales[self.density_unit]

        for law in mklaws(waves, types, modes, freqs / time_scale , values * dist_scale / time_scale, dvalues=None):
            key = "%s%s%d" % (law.wave.upper(), law.type.upper(), law.mode)
            coll = self.dispcoll[key]
            if self.period:
                coll['segments'].append(np.column_stack((1. / law.freq, law.value)))
            else:
                coll['segments'].append(np.column_stack((law.freq, law.value)))

            coll['colorvalues'].append(colorvalue)

    def showdepthcoll(self, vmin, vmax, cmap, **kwargs):
        for key, ax in self.axdepth.items():
            if ax is None :
                # compact
                continue
            segments = np.asarray(self.deptcoll[key]['segments'])
            if not len(segments):
                continue
            colorvalues = np.asarray(self.deptcoll[key]['colorvalues'])
            lc = LineCollection(segments, array=colorvalues, norm=plt.Normalize(vmin, vmax), cmap=cmap, **kwargs)
            ax.add_collection(lc)

    def showdispcoll(self, vmin, vmax, cmap, **kwargs):
        for key, ax in self.axdisp.items():
            coll = self.dispcoll[key]

            segments = coll['segments']
            colorvalues = coll['colorvalues']
            lc = LineCollection(segments, array=colorvalues, norm=plt.Normalize(vmin, vmax), cmap=cmap, **kwargs)
            ax.add_collection(lc)


# _____________________________________________
class DepthDispDisplayCompact(DepthDispDisplay):
    """display only vs and the dispersion curves"""
    def __init__(self, fig=None, targetfile=None):

        if fig is None:
            self.fig = plt.figure(figsize=(5, 6))#figsize=(18, 10))
            self.fig.subplots_adjust(wspace=0.05)
        else:
            fig.clf()
            self.fig = fig

        self.axdepth = {}
        self.axdepth['VS'] = self.fig.add_subplot(1, 2, 1, title="$V_S\,(%s/%s)$" % (self.dist_unit, self.time_unit), ylabel="depth (%s)" % self.dist_unit)
        self.axdepth['VP'] = self.axdepth['PR'] = self.axdepth['RH'] = None

        if targetfile is None:
            self.axdisp = {}
            self.axdisp['RU0'] = self.fig.add_subplot(3, 2, 2, title="RU0", ylabel="grpvel (%s/%s)" % (self.dist_unit, self.time_unit))
            self.axdisp['RU1'] = self.fig.add_subplot(3, 2, 4, title="RU1", ylabel="grpvel (%s/%s)" % (self.dist_unit, self.time_unit), sharex=self.axdisp['RU0'], sharey=self.axdisp['RU0'])

            for _, ax in self.axdisp.items():
                ax.loglog([1.],[1.]) #fuck matplotlib v2
                ax.yaxis.set_label_position("right")
                ax.yaxis.tick_right()
        else:
            # adjust the suplots according to the target data
            s = surf96reader(targetfile)
            ncurves = len(list(s.wtm()))
            ndiv = max([3, ncurves])

            self.axdisp = {}
            if ncurves <= 4:
                share = None
                wtms = list(s.wtm())
                for n, (w, t, m)  in enumerate(wtms):
                    ax = self.fig.add_subplot(ndiv, 2, (n+1)*2,
                                              title="%s%s%d" % (w, t, m),
                                              sharex=share, sharey=share,
                                              #ylabel="velocity (km/s)")
                                              ylabel="%s (%s/%s)" % (["grpvel", "phsvel"][int(t=="C")],  self.dist_unit, self.time_unit))
                    ax.yaxis.set_label_position("right")
                    ax.yaxis.tick_right()

                    self.axdisp["%s%s%d" % (w.upper(), t.upper(), m)] = share = ax
                    if n < len(wtms)-1:
                        plt.setp(ax.get_xticklabels(), visible = False)

            else:
                # all curves on same graph
                ax = self.fig.add_subplot(3, 2, 2,
                                          ylabel="Velocity (%s/%s)" % (self.dist_unit, self.time_unit))
                ax.yaxis.set_label_position("right")
                ax.yaxis.tick_right()
                for n, (w, t, m)  in enumerate(s.wtm()):
                    self.axdisp["%s%s%d" % (w.upper(), t.upper(), m)] = ax

            if self.period:
                ax.set_xlabel('$Period\,(%s)$' % self.time_unit)
            else:
                ax.set_xlabel('$Frequency\,(%s^{-1})$' % self.time_unit)

        pos = ax.get_position()

        self.cax = self.fig.add_axes((pos.x0, np.max([0.12, pos.y0 - 0.12]), pos.width, 0.01))
        # plt.setp(self.axdepth['VS'].get_yticklabels(), visible=False)
        self.axconv = None #Not needed here

        if not self.axdepth['VS'].yaxis_inverted():
            self.axdepth['VS'].invert_yaxis()

        self.clear_collections()

    def scale_and_plot_depthmodel1D(self, dm1d, which, **kwargs):
        if which.lower() in ["vp", "pr", "rh"]:
            return
        DepthDispDisplay.scale_and_plot_depthmodel1D(
            self, dm1d, which=which, **kwargs)

    def plotmodel(self, *args, **kwargs):
        kwargs.setdefault('showpr', False)
        kwargs.setdefault('showvp', False)
        kwargs.setdefault('showrh', False)

        DepthDispDisplay.plotmodel(self, *args, **kwargs)

    # _____________________________________________
    def addmodel(self, *args, **kwargs):
        kwargs.setdefault('showpr', False)
        kwargs.setdefault('showvp', False)
        kwargs.setdefault('showrh', False)

        DepthDispDisplay.addmodel(self, *args, **kwargs)


class DepthDispDisplaySI(DepthDispDisplay):
    dist_unit = "m"
    time_unit = "s"
    density_unit = "kg/m3"  # "g/cm3"
    period = False
    
    def tick(self):
        for _, ax in self.axdisp.items():
            logtick(ax, "xy")
            plt.setp(ax.get_yticklabels(),
                     rotation=-90,
                     ha="left", va="center",
                     fontsize=4)
            ax.tick_params(axis='y', which='major', pad=0)

        for _, ax in self.axdepth.items():
            if ax is not None:
                Ntick(ax, 4, "x")


class DepthDispDisplayCompactSI(DepthDispDisplayCompact):
    dist_unit = "m"
    time_unit = "s"
    density_unit = "kg/m3"  # "g/cm3"
    period = False
    
    def tick(self):
        for _, ax in self.axdisp.items():
            logtick(ax, "xy")
            plt.setp(ax.get_yticklabels(),
                     rotation=-90,
                     ha="left", va="center",
                     fontsize=4)
            ax.tick_params(axis='y', which='major', pad=0)



        for _, ax in self.axdepth.items():
            if ax is not None:
                Ntick(ax, 4, "x")

