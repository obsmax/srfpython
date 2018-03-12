import os, glob
import numpy as np
from srfpython.standalone.multipro8 import Job, MapAsync
from srfpython.standalone.display import values2colors, makecolorbar, legendtext, chftsz
from srfpython.standalone import cmaps #used in function eval
from srfpython.standalone.asciifile import AsciiFile
from srfpython.Herrmann.Herrmann import groupbywtm, igroupbywtm, dispersion
from srfpython.depthdisp.depthmodels import depthmodel_from_mod96, depthmodel_from_arrays, depthmodel
from srfpython.depthdisp.dispcurves import freqspace
from srfpython.depthdisp.depthpdfs import dmstats1
from srfpython.depthdisp.disppdfs import dispstats
from srfpython.depthdisp.depthdispdisplay import DepthDispDisplay, plt, showme
from srfpython.HerrMet.files import load_paramfile, RunFile
from srfpython.HerrMet.datacoders import makedatacoder, Datacoder_log
from srfpython.HerrMet.theory import overdisp
from srfpython.HerrMet.parameterizers import Parameterizer_mZVSPRRH, Parameterizer_mZVSVPRH, \
    Parameterizer_mZVSPRzRHz, Parameterizer_mZVSPRzRHvp


# ------------------------------ defaults
default_rootnames = "_HerrMet_*"
default_top_llkmin = -1000.
default_top_limit = 100
default_top_step = 1
default_pdf_llkmin = 0
default_pdf_limit = 0
default_pdf_step = 1
default_cmap = "gray" # plt.cm.jet# plt.cm.gray #


# ------------------------------ autorized_keys
authorized_keys = ["-top", "-overdisp", "-pdf", "-png", "-m96", "-cmap", "-inline", "-ritt"]

# ------------------------------ help messages
short_help = "--display    display target, parameterization, solutions"

long_help = """\
--display   s [s...] display param, target, and run outputs for the required rootnames, default {default_rootnames}
                     (use "." to see the parameterzation template ./_HerrMet.param from option --param)
    -top    [f i i]  show the best models on the figure, see --extract for arguments
                     default {default_top_llkmin}, {default_top_limit}, {default_top_step}
    -overdisp        recompute dispersion curves of the best models selected with higher resolution
    -pdf    [f i i]  compute and show the statistics for the selected models, see --extract for arguments
                     default {default_pdf_llkmin}, {default_pdf_limit}, {default_pdf_step} 
                     use --extract to save pdf outputs (median 16% and 84% percentiles at mod96 format)
    -png             save figure as pngfile instead of displaying it on screen
    -m96    s [s...] append depth model(s) to the plot from mod96 file(s)
    -cmap            colormap, default {default_cmap}
    -inline          do not pause (use in jupyter notebooks)
    """.format(default_rootnames=default_rootnames,
               default_top_llkmin=default_top_llkmin,
               default_top_limit=default_top_limit,
               default_top_step=default_top_step,
               default_pdf_llkmin=default_pdf_llkmin,
               default_pdf_limit=default_pdf_limit,
               default_pdf_step=default_pdf_step,
               default_cmap=default_cmap)

# ------------------------------ example usage
example = """\
## DISPLAY
# display the best 10 models,
# recompute the disp curves with higher resolution
# compute median and percentiles over the whole population of models
# save as png file (non display backend) 

HerrMet --display \\
            -top 0.0 10 1 \\
            -overdisp \\
            -pdf \\
            -png 
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
    if os.path.exists(targetfile):
        rd = DepthDispDisplay(targetfile=targetfile)
        d = makedatacoder(targetfile, which=Datacoder_log)  # datacoder based on observations
        dobs, _ = d.target()
    else:
        print "no target file found in %s" % rootname
        rd = DepthDispDisplay()

    # ------ Display run results if exist
    if os.path.exists(runfile) and ("-top" in argv.keys() or "-pdf" in argv.keys()):

        with RunFile(runfile, verbose=verbose) as rundb:
            s = rundb.select('select MODELID from MODELS limit 1')
            if s is not None:
                # --- display best models
                if "-top" in argv.keys():

                    assert argv["-top"] == [] or len(argv["-top"]) == 3  # unexpected argument number
                    if argv["-top"] == []:
                        top_llkmin, top_limit, top_step = default_top_llkmin, default_top_limit, default_top_step
                    elif len(argv['-top']) == 3:
                        top_llkmin, top_limit, top_step = argv['-top']

                    print "top : llkmin %f, limit %d, step %d" % (top_llkmin, top_limit, top_step),
                    chainids, weights, llks, ms, ds = rundb.getzip(llkmin=top_llkmin,
                                                                   limit=top_limit,
                                                                   step=top_step,
                                                                   algo="METROPOLIS")
                    vmin, vmax = llks.min(), llks.max()
                    colors = values2colors(llks, vmin=vmin, vmax=vmax, cmap=argv['-cmap'])

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
                        for clr, (mms, dds) in zip(colors[::-1],
                                                   overdisp(ms[::-1],
                                                            overwaves, overtypes, overmodes, overfreqs,
                                                            verbose=verbose, **mapkwargs)):
                            rd.plotmodel(color=clr, alpha=1.0, linewidth=3, *mms)
                            try:
                                rd.plotdisp(color=clr, alpha=1.0, linewidth=3, *dds)
                            except KeyboardInterrupt:
                                raise
                            except Exception as e:
                                print "Error : could not plot dispersion curve (%s)" % str(e)

                        cb = makecolorbar(vmin=vmin, vmax=vmax, cmap=argv['-cmap'])
                        pos = rd.axdisp[-1].get_position()
                        cax = rd.fig.add_axes((pos.x0, 0.12, pos.width, 0.01))
                        rd.fig.colorbar(cb, cax=cax, label="log likelyhood", orientation="horizontal")
                        cax.set_xticklabels(cax.get_xticklabels(), rotation=90., horizontalalignment="center")

                    else:
                        "display the dispersion curves as stored in the database"
                        for i in range(len(llks))[::-1]:
                            rd.plotmodel(color=colors[i], alpha=1.0, linewidth=3, *ms[i])
                            rd.plotdisp(color=colors[i], alpha=1.0, linewidth=3, *ds[i])

                        cb = makecolorbar(vmin=vmin, vmax=vmax, cmap=argv['-cmap'])
                        pos = rd.axdisp[-1].get_position()
                        cax = rd.fig.add_axes((pos.x0, 0.12, pos.width, 0.01))
                        rd.fig.colorbar(cb, cax=cax, label="log likelyhood", orientation="horizontal")
                        cax.set_xticklabels(cax.get_xticklabels(), rotation=90., horizontalalignment="center")

                # ---- display posterior pdf
                if "-pdf" in argv.keys():

                    assert argv["-pdf"] == [] or len(argv["-pdf"]) == 3  # unexpected argument number
                    if argv["-pdf"] == []:
                        pdf_llkmin, pdf_limit, pdf_step = default_pdf_llkmin, default_pdf_limit, default_pdf_step
                    elif len(argv['-pdf']) == 3:
                        pdf_llkmin, pdf_limit, pdf_step = argv['-pdf']

                    print "pdf : llkmin %f, limit %d, step %d" % (pdf_llkmin, pdf_limit, pdf_step),
                    chainids, weights, llks, ms, ds = rundb.getzip(llkmin=pdf_llkmin,
                                                                   limit=pdf_limit,
                                                                   step=pdf_step,
                                                                   algo="METROPOLIS")

                    dms = [depthmodel_from_arrays(ztop, vp, vs, rh) for ztop, vp, vs, rh in ms]

                    # # prepare the forward projection of the depth models
                    # fwdwaves, fwdtypes, fwdmodes, fwdfreqs = \
                    #     [np.array(_) for _ in d.waves, d.types, d.modes, d.freqs]
                    # if "-overdisp" in argv.keys():
                    #     overwaves, overtypes, overmodes, _, _ = zip(
                    #         *list(groupbywtm(fwdwaves, fwdtypes, fwdmodes, fwdfreqs, np.arange(len(fwdfreqs)), None, True)))
                    #     overfreqs = np.array([freqspace(0.6 * min(fwdfreqs), 1.4 * max(fwdfreqs), 100, "plog") for _ in
                    #                  xrange(len(overwaves))])
                    #     fwdwaves, fwdtypes, fwdmodes, fwdfreqs = \
                    #         igroupbywtm(overwaves, overtypes, overmodes, overfreqs)

                    # display the depth models and their projections
                    for p, (vppc, vspc, rhpc, prpc) in \
                            dmstats1(dms,
                                     percentiles=[0.01, 0.16, 0.5, 0.84, 0.99],
                                     Ndepth=100,
                                     Nvalue=100,
                                     weights=weights,
                                     **mapkwargs):
                        try:
                            l = 3 if p == 0.5 else 1
                            vppc.show(rd.axvp, color="b", linewidth=l)
                            vspc.show(rd.axvs, color="b", linewidth=l)
                            rhpc.show(rd.axrh, color="b", linewidth=l)
                            prpc.show(rd.axpr, color="b", linewidth=l)

                            # # project these solutions to the dataspace
                            # fwdvalues = dispersion(vppc.z, vppc.values, vspc.values, rhpc.values,
                            #                        fwdwaves, fwdtypes, fwdmodes, fwdfreqs)
                            #
                            # rd.plotdisp(fwdwaves, fwdtypes, fwdmodes, fwdfreqs, fwdvalues,
                            #              dvalues=None, color="b", alpha=1., linewidth=l)

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
                                        dvalues=None, color="b", alpha=1., linewidth=l)

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
        vplow.show(rd.axvp, alpha=1.0, color="r", marker="o--", linewidth=3)
        vphgh.show(rd.axvp, alpha=1.0, color="r", marker="o--", linewidth=3)
        vslow.show(rd.axvs, alpha=1.0, color="r", marker="o--", linewidth=3)
        vshgh.show(rd.axvs, alpha=1.0, color="r", marker="o--", linewidth=3)
        rhlow.show(rd.axrh, alpha=1.0, color="r", marker="o--", linewidth=3)
        rhhgh.show(rd.axrh, alpha=1.0, color="r", marker="o--", linewidth=3)
        prlow.show(rd.axpr, alpha=1.0, color="r", marker="o--", linewidth=3)
        prhgh.show(rd.axpr, alpha=1.0, color="r", marker="o--", linewidth=3)
        zmax = 1.1 * p.inv(p.MINF)[0][-1]

        if isinstance(p, Parameterizer_mZVSPRzRHvp):
            rd.axpr.plot(
                p.PRz(np.linspace(0., zmax, 100)),
                np.linspace(0., zmax, 100), "r--", linewidth=3)
            legendtext(rd.axpr, p.PRzName, loc=4)
            legendtext(rd.axrh, p.RHvpName, loc=4)
        elif isinstance(p, Parameterizer_mZVSPRzRHz):
            rd.axpr.plot(
                p.PRz(np.linspace(0., zmax, 100)),
                np.linspace(0., zmax, 100), "r--", linewidth=3)
            rd.axrh.plot(
                p.RHz(np.linspace(0., zmax, 100)),
                np.linspace(0., zmax, 100), "r--", linewidth=3)
            legendtext(rd.axpr, p.PRzName, loc=4)
            legendtext(rd.axrh, p.RHzName, loc=4)

        rd.set_zlim(np.array([0, zmax]))
    else:
        print "call option --param to see prior depth boundaries"

    # --------------------
    if "-m96" in argv.keys():  # plot user data on top
        for m96 in argv['-m96']:
            try:
                dm = depthmodel_from_mod96(m96)
                dm.vp.show(rd.axvp, "m", linewidth=3, label=m96)
                dm.vs.show(rd.axvs, "m", linewidth=3)
                dm.rh.show(rd.axrh, "m", linewidth=3)
                dm.pr().show(rd.axpr, "m", linewidth=3)
            except KeyboardInterrupt:
                raise
            except :#Exception as e:
                print 'could not read or display %s (reason : %s)' % (m96, str(e))
            rd.axvp.legend(loc=3)
    if "-ritt" in argv.keys():
        a = AsciiFile('/mnt/labex2/home/max/data/boreholes/GRT1/GRT1.logsonic')
        rd.axvs.plot(a.data['VS'], a.data['TVD']/1000., "m", alpha=0.5)
        rd.axvp.plot(a.data['VP'], a.data['TVD'] / 1000., "m", alpha=0.5)
        rd.axpr.plot(a.data['VP']/a.data['VS'], a.data['TVD'] / 1000., "m", alpha=0.5)
    # --------------------
    if os.path.exists(targetfile):
        # plot data on top
        rd.plotdisp(d.waves, d.types, d.modes, d.freqs, d.inv(dobs),
                    dvalues=d.dvalues,
                    alpha=0.8,
                    color="g",
                    linewidth=1)

        if "-overdisp" in argv.keys():
            rd.set_vlim((0.5 * d.values.min(), 1.5 * d.values.max()))
            rd.set_plim((0.8 / overfreqs.max(), 1.2 / overfreqs.min()))
        else:
            rd.set_vlim((0.8 * d.values.min(), 1.2 * d.values.max()))
            rd.set_plim((0.8 / d.freqs.max(), 1.2 / d.freqs.min()))
    rd.tick()
    rd.grid()
    rd.fig.suptitle(rootname.split('_HerrMet_')[-1])
    chftsz(rd.fig, 10)
    if "-png" in argv.keys():
        if verbose:
            print "writing %s" % pngfile
        rd.fig.savefig(pngfile)
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