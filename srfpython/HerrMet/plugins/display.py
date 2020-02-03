import os, glob
import numpy as np
from srfpython.standalone.multipro8 import Job, MapAsync
from srfpython.standalone.display import values2colors, legendtext, chftsz
from srfpython.standalone import cmaps #used in function eval
from srfpython.standalone.asciifile import AsciiFile
from srfpython.Herrmann.Herrmann import groupbywtm, igroupbywtm, dispersion
from srfpython.depthdisp.depthmodels import depthmodel_from_mod96, depthmodel_from_arrays, depthmodel
from srfpython.depthdisp.dispcurves import freqspace
from srfpython.depthdisp.depthpdfs import dmstats1
from srfpython.depthdisp.disppdfs import dispstats
from srfpython.depthdisp.depthdispdisplay import DepthDispDisplay, DepthDispDisplayCompact, plt, showme
from srfpython.HerrMet.files import load_paramfile, RunFile
from srfpython.HerrMet.datacoders import makedatacoder, Datacoder_log
from srfpython.HerrMet.theory import overdisp
from srfpython.HerrMet.parameterizers import Parameterizer_mZVSPRRH, Parameterizer_mZVSVPRH, \
    Parameterizer_mZVSPRzRHz, Parameterizer_mZVSPRzRHvp


# ------------------------------ defaults
default_rootnames = "_HerrMet_*"
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
authorized_keys = ["-plot", "-overdisp", "-pdf", "-png", "-m96", "-cmap", "-compact", "-ftsz", "-inline", "-ritt"]

# ------------------------------ help messages
short_help = "--display    display target, parameterization, solutions"

long_help = """\
--display   s [s...] display param, target, and run outputs for the required rootnames, default {default_rootnames}
                     (use "." to see the parameterzation template ./_HerrMet.param from option --param)
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
    -png   [i]       save figure as pngfile instead of displaying it on screen, requires dpi, default {default_dpi} 
    -m96    s [s...] append depth model(s) to the plot from mod96 file(s)
    -cmap            colormap, default {default_cmap}
    -compact         display only vs and the dispersion curves, default False
    -ftsz  i         set font size, default {default_fontsize}
    -inline          do not pause (use in jupyter notebooks)
    """.format(default_rootnames=default_rootnames,
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
def _display_function(rootname, argv, verbose, mapkwargs):
    """private"""
    targetfile = "%s/_HerrMet.target" % rootname
    paramfile = "%s/_HerrMet.param" % rootname
    runfile = '%s/_HerrMet.run' % rootname
    pngfile = '%s/_HerrMet.png' % rootname
    #HerrLininitfile = '%s/_HerrLin.init' % rootname

    # ------ Initiate the displayer using the target data if exists
    if "-compact" in argv.keys(): # compact mode
        which_displayer = DepthDispDisplayCompact
    else:
        which_displayer = DepthDispDisplay

    if os.path.exists(targetfile):
        rd = which_displayer(targetfile=targetfile)
        d = makedatacoder(targetfile, which=Datacoder_log)  # datacoder based on observations
        dobs, _ = d.target()
    else:
        print "no target file found in %s" % rootname
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

                    print "plot : %s, limit %d, llkmin %f, step %d" % (plot_mode, plot_limit, plot_llkmin, plot_step),
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
                                     xrange(len(overwaves))]
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
                                print "Error : could not plot dispersion curve (%s)" % str(e)

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
                    print rd.cax.get_position()
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

                    print "pdf : %s, limit %d, llkmin %f, step %d" % (pdf_mode, pdf_limit, pdf_llkmin, pdf_step),
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
                            for what, where in zip([vppc, vspc, rhpc, prpc], [rd.axdepth['VP'], rd.axdepth['VS'], rd.axdepth['RH'], rd.axdepth['PR']]):
                                if where is not None:
                                    what.show(where, color=clr, linewidth=l, alpha=alp)

                        except KeyboardInterrupt:
                            raise
                        except Exception as e:
                            print "Error", str(e)

                    # display the disp pdf
                    for p, (wpc, tpc, mpc, fpc, vpc) in \
                            dispstats(ds,
                                      percentiles=[0.01, 0.16, 0.5, 0.84, 0.99],
                                      Ndisp=100,
                                      weights=weights,
                                      **mapkwargs):
                        try:
                            l = 3 if p == 0.5 else 1
                            rd.plotdisp(wpc, tpc, mpc, fpc, vpc,
                                        dvalues=None, color=clr, alpha=alp, linewidth=l)

                        except KeyboardInterrupt:
                            raise
                        except Exception as e:
                            print "Error", str(e)

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
        else:
            raise Exception('')

        #
        vplow, vphgh, vslow, vshgh, rhlow, rhhgh, prlow, prhgh = p.boundaries()

        for what, where in zip(\
                [vplow, vphgh, vslow, vshgh, rhlow, rhhgh, prlow, prhgh],
                [rd.axdepth['VP'], rd.axdepth['VP'], rd.axdepth['VS'], rd.axdepth['VS'], rd.axdepth['RH'], rd.axdepth['RH'], rd.axdepth['PR'], rd.axdepth['PR']]):
            if where is not None:
                what.show(where, alpha=1.0, color="k", marker="o--", linewidth=1, markersize=3)
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
        print "call option --param to see prior depth boundaries"

    # --------------------
    if "-m96" in argv.keys():  # plot user data on top
        for m96 in argv['-m96']:
            try:
                dm = depthmodel_from_mod96(m96)
                dm.vp.show(rd.axdepth['VP'], "m", linewidth=3, label=m96)
                dm.vs.show(rd.axdepth['VS'], "m", linewidth=3)
                dm.rh.show(rd.axdepth['RH'], "m", linewidth=3)
                dm.pr().show(rd.axdepth['PR'], "m", linewidth=3)
            except KeyboardInterrupt:
                raise
            except :#Exception as e:
                print 'could not read or display %s (reason : %s)' % (m96, str(e))
            rd.axdepth['VP'].legend(loc=3)
    if "-ritt" in argv.keys():
        a = AsciiFile('/mnt/labex2/home/max/data/boreholes/GRT1/GRT1.logsonic')

        for what, where in zip([a.data['VS'], a.data['VP'], a.data['VP']/a.data['VS']], [rd.axdepth['VS'], rd.axdepth['VP'], rd.axdepth['PR']]):
            if where is not None:
                where.plot(what, a.data['TVD']/1000., "m", alpha=0.5)

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
    rd.fig.suptitle(rootname.split('_HerrMet_')[-1])
    if "-ftsz" in argv.keys():
        chftsz(rd.fig, argv["-ftsz"][0])
    else:
        chftsz(rd.fig, default_fontsize)
    if "-png" in argv.keys():
        dpi = argv['-png'][0] if len(argv['-png']) else default_dpi
        if verbose:
            print "writing %s" % pngfile
        rd.fig.savefig(pngfile, dpi=dpi)
    elif "-inline" in argv.keys():
        plt.show()
    else:
        showme()
    plt.close(rd.fig)


# ------------------------------
def display(argv, verbose, mapkwargs):

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

    # ----------- special case, just show the parameterization file from --param : ./_HerrMet.param
    if len(rootnames) == 1 and rootnames[0] == '.':
        _display_function(".", argv=argv, verbose=verbose, mapkwargs=mapkwargs)

    # ----------- general case
    else:
        for rootname in rootnames:
            if not os.path.isdir(rootname):
                raise Exception('%s does not exist' % rootname)
            elif not rootname.startswith('_HerrMet_'):
                raise Exception('%s does not starts with _HerrMet_' % rootname)

        if "-png" not in argv.keys():
            # display mode, cannot parallelize
            for rootname in rootnames:
                _display_function(rootname, argv=argv, verbose=verbose, mapkwargs=mapkwargs)
        else:
            def gen():
                for rootname in rootnames:
                    yield Job(rootname, argv, verbose=verbose, mapkwargs=mapkwargs)

            with MapAsync(_display_function, gen(), **mapkwargs) as ma:
                for _ in ma:
                    pass