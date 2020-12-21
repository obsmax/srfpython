from __future__ import print_function

import os, glob
import numpy as np
from srfpython.standalone.multipro8 import Job, MapAsync
from srfpython.standalone.display import values2colors, legendtext, chftsz
from srfpython.standalone import cmaps #used in function eval
from srfpython.standalone.asciifile import AsciiFile
from srfpython.Herrmann.Herrmann import groupbywtm, igroupbywtm
from srfpython.depthdisp.depthmodels import depthmodel_from_mod96, depthmodel_from_arrays, depthmodel
from srfpython.depthdisp.dispcurves import freqspace
from srfpython.depthdisp.depthpdfs import dmstats1
from srfpython.depthdisp.disppdfs import dispstats
from srfpython.depthdisp.depthdispdisplay import DepthDispDisplay, DepthDispDisplayCompact, plt, showme
from srfpython.HerrMet.runfile import RunFile
from srfpython.HerrMet.paramfile import load_paramfile
from srfpython.HerrMet.datacoders import makedatacoder, Datacoder_log
from srfpython.HerrMet.overdisp import overdisp
from srfpython.HerrMet.parameterizers import Parameterizer_mZVSPRRH, Parameterizer_mZVSVPRH, \
    Parameterizer_mZVSPRzRHz, Parameterizer_mZVSPRzRHvp, Parameterizer_mZVSVPvsRHvp
from srfpython.HerrMet.files import \
    DEFAULTROOTNAMES, HERRMETTARGETFILE, HERRMETRUNFILE, HERRMETPARAMFILE, \
    HERRMETPARAMFILELOCAL, HERRMETDISPLAYFILE, rootname_to_nodename

# ------------------------------ defaults
default_rootnames = DEFAULTROOTNAMES

default_plot_mode = "last"
default_plot_limit = 100
default_plot_llkmin = 0.
default_plot_step = 1
default_pdf_mode = "last"
default_pdf_limit = 0
default_pdf_llkmin = 0.
default_pdf_step = 1
default_cmap = "viridis"  # plt.cm.jet# plt.cm.gray #
default_fontsize = 10
default_dpi = 100


# ------------------------------ autorized_keys
authorized_keys = ["-plot", "-overdisp", "-pdf", "-png", "-svg", "-m96", "-cmap", "-compact", "-ftsz", "-inline", "-h", "-help"]

# ------------------------------ help messages
short_help = "--display    display target, parameterization, solutions"

long_help = """\
--display   s [s...] display param, target, and run outputs for the required rootnames, default {default_rootnames}
                     (use "." to see the parameterzation template {herrmetparamfilelocal} from option --param)
    -plot  [s i f i] show the best models on the figure, arguments are :  
                     first argument = selection mode, last or best
                     second argument = highest model number to include (>=0, 0 means all)  
                     third argument = lowest log likelihood value to include (<=0.0, 0.0 means all)
                     fourth argument = include only one model over "step" (>=1)
                     default {default_plot_mode} {default_plot_limit} {default_plot_llkmin} {default_plot_step}             
    -overdisp        recompute dispersion curves of the best models selected with higher resolution
    -pdf   [s i f i] compute and show the statistics for the selected models, see -plot for arguments
                     default {default_pdf_mode} {default_pdf_limit} {default_pdf_llkmin} {default_pdf_step} 
                     use --extract to save pdf outputs
    -png/-svg [i]    save figure instead of displaying it on screen, requires dpi, default {default_dpi} 
    -m96    s [s...] append depth model(s) to the plot from mod96 file(s)
    -cmap            colormap, default {default_cmap}
    -compact         display only vs and the dispersion curves, default False
    -ftsz  i         set font size, default {default_fontsize}
    -inline          do not pause (use in jupyter notebooks)
    -h, -help        display the help message for this plugin     
    """.format(
               herrmetparamfilelocal=HERRMETPARAMFILELOCAL,
               default_rootnames=default_rootnames,
               default_plot_mode=default_plot_mode,
               default_plot_limit=default_plot_limit,
               default_plot_llkmin=default_plot_llkmin,
               default_plot_step=default_plot_step,
               default_pdf_mode=default_pdf_mode,
               default_pdf_limit=default_pdf_limit,
               default_pdf_llkmin=default_pdf_llkmin,
               default_pdf_step=default_pdf_step,
               default_cmap=default_cmap,
               default_fontsize=default_fontsize,
               default_dpi=default_dpi)

# ------------------------------ example usage
example = """\
## DISPLAY
# display the best 10 models,
# recompute the disp curves with higher resolution
# compute median and percentiles over the whole population of models
# save as png file (non display backend) 

HerrMet --display \\
            -plot last 10 0.0 1 \\
            -overdisp \\
            -pdf \\
            -png 300
"""


# -------------------------------------
def _display_function(rootname, argv, verbose, mapkwargs, fig=None, return_fig=False):
    """private"""
    targetfile = HERRMETTARGETFILE.format(rootname=rootname)
    paramfile = HERRMETPARAMFILE.format(rootname=rootname)
    runfile = HERRMETRUNFILE.format(rootname=rootname)


    # ------ Initiate the displayer using the target data if exists
    if "-compact" in argv.keys():  # compact mode
        which_displayer = DepthDispDisplayCompact
        displayfile = HERRMETDISPLAYFILE.format(rootname=rootname, options="_compact")
    else:
        which_displayer = DepthDispDisplay
        displayfile = HERRMETDISPLAYFILE.format(rootname=rootname, options="")

    if os.path.exists(targetfile):
        rd = which_displayer(targetfile=targetfile, fig=fig)
        d = makedatacoder(targetfile, which=Datacoder_log)  # datacoder based on observations
        dobs, _ = d.target()
    else:
        print("no target file found in %s" % rootname)
        rd = which_displayer()

    # ------ Display run results if exist
    if os.path.exists(runfile) and ("-plot" in argv.keys() or "-pdf" in argv.keys()):

        with RunFile(runfile, verbose=verbose) as rundb:
            s = rundb.select('select MODELID from MODELS limit 1')
            if s is not None:
                # --- display best models
                if "-plot" in argv.keys():

                    assert argv["-plot"] == [] or len(argv["-plot"]) == 4  # unexpected argument number
                    if argv["-plot"] == []:
                        plot_mode, plot_limit, plot_llkmin, plot_step = \
                            default_plot_mode, default_plot_limit, \
                            default_plot_llkmin, default_plot_step
                    elif len(argv['-plot']) == 4:
                        plot_mode, plot_limit, plot_llkmin, plot_step = argv['-plot']
                    else:
                        raise Exception()

                    print("plot : %s, limit %d, llkmin %f, step %d" %
                        (plot_mode, plot_limit, plot_llkmin, plot_step), end=' ')
                    if plot_mode == "best":
                        chainids, weights, llks, ms, ds = \
                            rundb.getzip(limit=plot_limit,
                                         llkmin=plot_llkmin,
                                         step=plot_step,
                                         algo="METROPOLIS")
                    elif plot_mode == "last":
                        chainids, weights, llks, ms, ds = \
                            rundb.getlastszip(limit=plot_limit,
                                              llkmin=plot_llkmin,
                                              step=plot_step,
                                              algo="METROPOLIS")
                    else:
                        raise Exception('unexpected plot mode %s' % plot_mode)

                    vmin, vmax = llks.min(), llks.max()
                    # colors = values2colors(llks, vmin=vmin, vmax=vmax, cmap=argv['-cmap'])

                    if "-overdisp" in argv.keys():
                        """note : recomputing dispersion with another frequency array might
                                  result in a completely different dispersion curve in case
                                  of root search failure """
                        waves, types, modes, freqs, _ = ds[0]
                        overwaves, overtypes, overmodes, _, _ = zip(
                            *list(groupbywtm(waves, types, modes, freqs, np.arange(len(freqs)), None, True)))
                        overfreqs = [freqspace(0.6 * min(freqs), 1.4 * max(freqs), 100, "plog") for _ in
                                     range(len(overwaves))]
                        overwaves, overtypes, overmodes, overfreqs = \
                            igroupbywtm(overwaves, overtypes, overmodes, overfreqs)
                        for llk, (mms, dds) in zip(llks[::-1],
                                                   overdisp(ms[::-1],
                                                            overwaves, overtypes, overmodes, overfreqs,
                                                            verbose=verbose, **mapkwargs)):
                            # rd.plotmodel(color=clr, alpha=1.0, linewidth=3, *mms)
                            rd.addmodel(colorvalue=llk, *mms)
                            try:
                                # rd.plotdisp(color=clr, alpha=1.0, linewidth=3, *dds)
                                rd.adddisp(colorvalue=llk, *dds)
                            except KeyboardInterrupt:
                                raise
                            except Exception as e:
                                print("Error : could not plot dispersion curve (%s)" % str(e))

                        # cb = makecolorbar(vmin=vmin, vmax=vmax, cmap=argv['-cmap'])
                        # pos = rd.axdisp[-1].get_position()
                        # cax = rd.fig.add_axes((pos.x0, 0.12, pos.width, 0.01))
                        # rd.fig.colorbar(cb, cax=cax, label="log likelihood", orientation="horizontal")


                    else:
                        "display the dispersion curves as stored in the database"
                        for i in range(len(llks))[::-1]:
                            # rd.plotmodel(color=colors[i], alpha=1.0, linewidth=3, *ms[i])
                            # rd.plotdisp(color=colors[i], alpha=1.0, linewidth=3, *ds[i])
                            rd.addmodel(colorvalue=llks[i], *ms[i])
                            rd.adddisp(colorvalue=llks[i], *ds[i])
                        # cb = makecolorbar(vmin=vmin, vmax=vmax, cmap=argv['-cmap'])
                        # pos = rd.axdisp[-1].get_position()
                        # cax = rd.fig.add_axes((pos.x0, 0.12, pos.width, 0.01))
                        # rd.fig.colorbar(cb, cax=cax, label="log likelihood", orientation="horizontal")
                        # cax.set_xticklabels(cax.get_xticklabels(), rotation=90., horizontalalignment="center")

                    rd.showdispcoll(vmin=vmin, vmax=vmax, cmap=argv['-cmap'], alpha=1.0, linewidth=3)
                    rd.showdepthcoll(vmin=vmin, vmax=vmax, cmap=argv['-cmap'], alpha=1.0, linewidth=3)
                    rd.colorbar(vmin=vmin, vmax=vmax, cmap=argv['-cmap'], label="log likelihood", orientation="horizontal")
                    # print rd.cax.get_position()
                    rd.cax.set_xticklabels(rd.cax.get_xticklabels(), rotation=90., horizontalalignment="center")

                # ---- display posterior pdf
                if "-pdf" in argv.keys():

                    assert argv["-pdf"] == [] or len(argv["-pdf"]) == 4  # unexpected argument number
                    if argv["-pdf"] == []:
                        pdf_mode, pdf_limit, pdf_llkmin, pdf_step = \
                            default_pdf_mode, default_pdf_limit, default_pdf_llkmin, default_pdf_step
                    elif len(argv['-pdf']) == 4:
                        pdf_mode, pdf_limit, pdf_llkmin, pdf_step = argv['-pdf']
                    else:
                        raise Exception()

                    print("pdf : %s, limit %d, llkmin %f, step %d" %
                          (pdf_mode, pdf_limit, pdf_llkmin, pdf_step), end=' ')
                    if pdf_mode == "best":
                        chainids, weights, llks, ms, ds = \
                            rundb.getzip(limit=pdf_limit,
                                         llkmin=pdf_llkmin,
                                         step=pdf_step,
                                         algo="METROPOLIS")
                    elif pdf_mode == "last":
                        chainids, weights, llks, ms, ds = \
                            rundb.getlastszip(limit=pdf_limit,
                                              llkmin=pdf_llkmin,
                                              step=pdf_step,
                                              algo="METROPOLIS")
                    else:
                        raise Exception('unexpected pdf mode %s' % pdf_mode)

                    dms = [depthmodel_from_arrays(ztop, vp, vs, rh) for ztop, vp, vs, rh in ms]

                    # display percentiles of model and data pdfs
                    clr = "b" if "-plot" not in argv.keys() else "k"
                    alp = 1.0 if "-plot" not in argv.keys() else 0.5

                    for p, (vppc, vspc, rhpc, prpc) in \
                            dmstats1(dms,
                                     percentiles=[0.01, 0.16, 0.5, 0.84, 0.99],
                                     Ndepth=100,
                                     Nvalue=100,
                                     weights=weights,
                                     **mapkwargs):
                        try:
                            l = 3 if p == 0.5 else 1
                            for what, where in zip(
                                    [vppc, vspc, rhpc, prpc],
                                    [rd.axdepth['VP'], rd.axdepth['VS'], rd.axdepth['RH'], rd.axdepth['PR']]):
                                if where is not None:
                                    what.show(where, color=clr, linewidth=l, alpha=alp)

                        except KeyboardInterrupt:
                            raise
                        except Exception as e:
                            print("Error", str(e))

                    # display the disp pdf
                    for p, (wpc, tpc, mpc, fpc, vpc) in \
                            dispstats(ds,
                                      percentiles=[0.01, 0.16, 0.5, 0.84, 0.99],
                                      Ndisp=200,
                                      weights=weights,
                                      **mapkwargs):
                        try:
                            l = 3 if p == 0.5 else 1
                            rd.plotdisp(wpc, tpc, mpc, fpc, vpc,
                                        dvalues=None, color=clr, alpha=alp, linewidth=l)

                        except KeyboardInterrupt:
                            raise
                        except Exception as e:
                            print("Error", str(e))

    # ------
    if os.path.exists(paramfile):
        p, _ = load_paramfile(paramfile)
        showvp, showvs, showrh, showpr = True, True, True, True
        if isinstance(p, Parameterizer_mZVSVPRH):
            showpr = False
        elif isinstance(p, Parameterizer_mZVSPRRH):
            showvp = False
        elif isinstance(p, Parameterizer_mZVSPRzRHvp):
            showvp = showpr = showrh = False
        elif isinstance(p, Parameterizer_mZVSPRzRHz):
            showvp = showpr = showrh = False
        elif isinstance(p, Parameterizer_mZVSVPvsRHvp):
            showvp = showpr = showrh = False
        else:
            raise NotImplementedError('type {} not implemented'.format(type(p)))

        #
        vplow, vphgh, vslow, vshgh, rhlow, rhhgh, prlow, prhgh = p.boundaries()

        for what, where in zip(
                [vplow, vphgh, vslow, vshgh, rhlow, rhhgh, prlow, prhgh],
                [rd.axdepth['VP'], rd.axdepth['VP'],
                 rd.axdepth['VS'], rd.axdepth['VS'],
                 rd.axdepth['RH'], rd.axdepth['RH'],
                 rd.axdepth['PR'], rd.axdepth['PR']]):

            if where is not None:
                what.show(
                    where,
                    alpha=1.0, color="k", marker="o--",
                    linewidth=1, markersize=3)

        vpmean, vsmean, prmean, rhmean = p.meanmodel()
        for what, where in zip(
                [vpmean, vsmean, prmean, rhmean],
                [rd.axdepth['VP'], rd.axdepth['VS'],
                 rd.axdepth['PR'], rd.axdepth['RH']]):

            if where is not None:
                what.show(
                    where,
                    alpha=1.0, color="r", marker="o--",
                    linewidth=1, markersize=3)

        zmax = 1.1 * p.inv(p.MINF)[0][-1]

        if isinstance(p, Parameterizer_mZVSPRzRHvp):
            rd.axdepth['PR'].plot(
                p.PRz(np.linspace(0., zmax, 100)),
                np.linspace(0., zmax, 100), "r--", linewidth=3)
            legendtext(rd.axdepth['PR'], p.PRzName, loc=4)
            legendtext(rd.axdepth['RH'], p.RHvpName, loc=4)
        elif isinstance(p, Parameterizer_mZVSPRzRHz):
            rd.axdepth['PR'].plot(
                p.PRz(np.linspace(0., zmax, 100)),
                np.linspace(0., zmax, 100), "r--", linewidth=3)
            rd.axdepth['RH'].plot(
                p.RHz(np.linspace(0., zmax, 100)),
                np.linspace(0., zmax, 100), "r--", linewidth=3)
            legendtext(rd.axdepth['PR'], p.PRzName, loc=4)
            legendtext(rd.axdepth['RH'], p.RHzName, loc=4)

        rd.set_zlim(np.array([0, zmax]))
    else:
        print("call option --param to see prior depth boundaries")

    # --------------------
    if "-m96" in argv.keys():  # plot user data on top
        for m96 in argv['-m96']:
            try:
                dm = depthmodel_from_mod96(m96)
                if "-compact" in argv.keys():  # compact mode 
                     dm.vs.show(rd.axdepth['VS'], "m", linewidth=3)               
                     rd.axdepth['VS'].legend(loc=3)
                     
                else:
                    dm.vp.show(rd.axdepth['VP'], "m", linewidth=3, label=m96)
                    dm.vs.show(rd.axdepth['VS'], "m", linewidth=3)
                    dm.rh.show(rd.axdepth['RH'], "m", linewidth=3)
                    dm.pr().show(rd.axdepth['PR'], "m", linewidth=3)
                    rd.axdepth['VP'].legend(loc=3)

            except KeyboardInterrupt:
                raise
            except Exception as e:
                print('could not read or display %s (reason : %s)' % (m96, str(e)))

    #     for what, where in zip([a.data['VS'], a.data['VP'], a.data['VP']/a.data['VS']], [rd.axdepth['VS'], rd.axdepth['VP'], rd.axdepth['PR']]):
    #         if where is not None:
    #             where.plot(what, a.data['TVD']/1000., "m", alpha=0.5)

    # --------------------
    if os.path.exists(targetfile):
        # plot data on top
        rd.plotdisp(d.waves, d.types, d.modes, d.freqs, d.inv(dobs),
                    dvalues=d.dvalues,
                    alpha=.5,
                    color="r",
                    linewidth=2)

        if "-overdisp" in argv.keys():
            rd.set_vlim((0.5 * d.values.min(), 1.5 * d.values.max()))
            rd.set_plim((0.8 / overfreqs.max(), 1.2 / overfreqs.min()))
        else:
            rd.set_vlim((0.8 * d.values.min(), 1.2 * d.values.max()))
            rd.set_plim((0.8 / d.freqs.max(), 1.2 / d.freqs.min()))
    rd.tick()
    rd.grid()
    rd.fig.suptitle(rootname_to_nodename(rootname))
    if "-ftsz" in argv.keys():
        chftsz(rd.fig, argv["-ftsz"][0])
    else:
        chftsz(rd.fig, default_fontsize)
    if "-png" in argv.keys() or "-svg" in argv.keys():

        if "-png" in argv.keys():
            dpi = argv['-png'][0] if len(argv['-png']) else default_dpi
            displayfile = displayfile.replace('.svg', ".png")
        elif "-svg" in argv.keys():
            dpi = argv['-svg'][0] if len(argv['-svg']) else default_dpi
            displayfile = displayfile.replace('.png', ".svg")
        if verbose:
            print("writing %s" % displayfile)
        rd.fig.savefig(displayfile, dpi=dpi)
    elif "-inline" in argv.keys():
        plt.show()
    else:
        showme()

    if return_fig:
        return rd.fig  # caution not for parallel apps
    else:
        plt.close(rd.fig)


# ------------------------------
def display(argv, verbose, mapkwargs):

    if '-h' in argv.keys() or "-help" in argv.keys():
        print(long_help)
        return

    for k in argv.keys():
        if k in ['main', "_keyorder"]:
            continue  # private keys

        if k not in authorized_keys:
            raise Exception('option %s is not recognized' % k)

    rootnames = argv['main']
    if rootnames == []:
        rootnames = glob.glob(default_rootnames)
    assert len(rootnames)

    # -------------------------------------
    if "-cmap" not in argv.keys():
        argv['-cmap'] = [default_cmap]

    try:
        argv['-cmap'] = plt.get_cmap(argv['-cmap'][0])
    except ValueError:
        try:
            argv['-cmap'] = eval("cmaps.%s()" % argv['-cmap'][0])
        except:
            raise Exception('could not find colormap %s neither in matplotlib neither in srfpython.standalone.utils.cmaps' % argv['-cmap'][0])

    # ----------- special case, just show the parameterization file from --param
    if len(rootnames) == 1 and rootnames[0] == '.':
        _display_function(".", argv=argv, verbose=verbose, mapkwargs=mapkwargs)

    # ----------- general case
    else:
        for rootname in rootnames:
            if not os.path.isdir(rootname):
                raise IOError('{} does not exist'.format(rootname))

        if "-png" not in argv.keys():
            # display mode, cannot parallelize
            fig = None
            for rootname in rootnames:
                fig = _display_function(
                    rootname, argv=argv, verbose=verbose, mapkwargs=mapkwargs, fig=fig, return_fig=True)
        else:
            def gen():
                for rootname in rootnames:
                    yield Job(rootname, argv, verbose=verbose, mapkwargs=mapkwargs, return_fig=False)

            with MapAsync(_display_function, gen(), **mapkwargs) as ma:
                for _ in ma:
                    pass
