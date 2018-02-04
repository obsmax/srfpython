#!/usr/bin/env python

#main imports
import sys, matplotlib
if "-png" in sys.argv[1:]: matplotlib.use('agg')
import os, imp, glob
import numpy as np

#tetedenoeud
from tetedenoeud.multipro.multipro8 import Job, MapAsync, MapSync
from tetedenoeud.utils.display import values2colors, makecolorbar, legendtext, chftsz
from tetedenoeud.utils import cmaps
from tetedenoeud.utils.cmaps import tej #to be evaluated
# from tetedenoeud.utils.asciifile import AsciiFile

#srfpython
from srfpython.Herrmann.Herrmann import groupbywtm, igroupbywtm, check_herrmann_codes
from srfpython.inversion.metropolis2 import LogGaussND, metropolis
from srfpython.depthdisp.depthmodels import depthmodel, depthmodel1D, depthmodel_from_mod96, depthmodel_from_arrays, depthspace
from srfpython.depthdisp.dispcurves import Claw, surf96reader, freqspace
from srfpython.depthdisp.depthpdfs import dmstats1
from srfpython.depthdisp.depthdispdisplay import DepthDispDisplay, plt, gcf, gca, showme, pause
from srfpython.utils import readargv, minmax, tostr

#local imports
from srfpython.HerrMet.files import write_default_paramfile, load_paramfile, RunFile
from srfpython.HerrMet.datacoders import makedatacoder, Datacoder_log
from srfpython.HerrMet.theory import Theory, overdisp
from srfpython.HerrMet.parameterizers import Parameterizer_mZVSPRRH, Parameterizer_mZVSVPRH, Parameterizer_mZVSPRzRHz, Parameterizer_mZVSPRzRHvp
#

check_herrmann_codes()
# -------------------------------------
version = "6.0"
default_verbose=True
default_rootnames = "_HerrMet_*"
default_nchain = 12
default_nkeep = 100
default_top_llkmin = -1000.
default_top_limit = 100
default_top_step = 1
default_pdf_llkmin = 0
default_pdf_limit = 0
default_pdf_step = 1
default_extract_llkmin = 0
default_extract_limit = 0
default_extract_step = 1
default_cmap = "gray" # plt.cm.jet# plt.cm.gray #
default_parameterization_list = ['mZVSPRRH', 'mZVSVPRH', 'mZVSPRzRHvp', 'mZVSPRzRHz']
default_parameterization = default_parameterization_list[0]
mapkwargs = {"Nworkers" : None, "Taskset" : None} #keyword arguments for every parallelized process
# -------------------------------------
autorizedkeys = \
    ["w", "taskset", "lowprio", "inline", "cmap",
     "version", "v",
     "help", "h",
     "example", "ex",
     "target", "resamp", "lunc", "unc", "ot",
     "param", "basedon", "t", "dvp", "dvs", "drh", "growing", "op",
     "send", "op",
     "run", "append", "nchain", "nkeep", "verbose",
     "extract",
     "display", "top", "overdisp", "pdf", "png", "m96",
     "test"]
# -------------------------------------
help = '''HerrMet V{version}
# ----------------------------------------------------
-w           i       set the number of virtual workers to use for all parallelized processes, default {default_Nworkers}
-taskset     s       change job affinity for all parallelized processes, default {default_Taskset}
-agg                 use agg backend (no display) if mentioned
-lowprio             run processes with low priority if mentioned
-inline              replace showme by plt.show (e.g. jupyter)
-verbose off         reduce verbosity
-cmap                colormap, default {default_cmap}
# ----------------------------------------------------
--version, -v        display version and quit
--help, -h           display this help message, and quit
--example, -ex       display an example of script, and quit
--target     s [s..] set the target dispersion curve(s) from surf96 file(s) (not modified)
                     for each target, I create a directory in . for temporary files
                     the data will be reproduced into a target file that can be customized manually
                     (to remove unwanted points, resample dispersion curves...) 
    -resamp  f f i s resample the dispersion curve in the target file, 
                     requires fmin(Hz), fmax(Hz), nfreq, fscale 
                     (flin=linear in freq domain, plin=linear in period or log=logarithmic scaling)
    -lunc    f       set constant uncertainties in log domain (uncertainty = value x lunc)
    -unc     f       set constant uncertainty in linear domain (uncertainty = unc)
    -ot              force overwriting _HerrMet.target if exists
--param      i f     generate a template parameter file to custom in .
                     need the number of layers and bottom depth in km
    -basedon s       build parametrization based on an existing mod96 file, require a filename, 
                     if not specified, I take fixed values to build the parameter file
    -t       s       parameterization type to use ({default_parameterization_list}), 
                     default {default_parameterization}
          mZVSPRRH = parameterize with 
                     - depth interface (mZ1 = first interface, mZ2, ...), 
                     - VS in each layer (VS0 = first layer, ...), 
                     - VP/VS in each layer (PR0, PR1, ...), 
                     - Density in each layer (RH0, RH1, ...)
          mZVSVPRH = parameterize with depth interface, VP in each layer, 
                     - depth interface (mZ1 = first interface, mZ2, ...), 
                     - VS in each layer (VS0 = first layer, ...), 
                     - VP in each layer (VP0, VP1, ...), 
                     - Density in each layer (RH0, RH1, ...)
       mZVSPRzRHvp = parameterize with depth interface, use fixed 
                     - depth interface (mZ1 = first interface, mZ2, ...), 
                     - VS in each layer (VS0 = first layer, ...), 
                     - use a fixed relation between VP/VS = f(z)
                     - use a fixed relation between RH = f(VP)
        mZVSPRzRHz = parameterize with depth interface, use fixed 
                     - depth interface (mZ1 = first interface, mZ2, ...), 
                     - VS in each layer (VS0 = first layer, ...), 
                     - use a fixed relation between VP/VS = f(z)
                     - use a fixed relation between RH = f(z) 
    -dvp     f f     add prior constraint on the vp offset between layers, 
                     requires the extremal values, km/s
    -dvs     f f     add prior constraint on the vs offset between layers, idem
    -drh     f f     add prior constraint on the density offset between layers, idem, g/cm3
    -dpr     f f     add prior constraint on the vp/vs offset between layers, idem, no unit
    -growing         shortcut for -dvp 0. 5. -dvs 0. 5. -drh 0. 5. -dpr -5. 0.
    -op              force overwriting ./_HerrMet.param if exists                     
--send       s [s..] send the custom parameterization file ./_HerrMet.param to the specified rootnames, 
                     default {default_rootnames}
    -op              force overwriting ./rootname/_HerrMet.param if exists
--run        s [s..] start inversion for the required rootnames, default {default_rootnames} 
    -append          set running mode to append (i.e. add new models to the exiting run file(s)) 
                     default is restart (overwrite the current run file(s) if any)
    -nchain  i       number of chains to use, default {default_nchain}
    -nkeep   i       number of models to keep per chain, default {default_nkeep}
    -w       i       see above, controls the max number of chains to run simultaneously
--extract    s [s..] rootnames for which to compute and extract pdf, default {default_rootnames}
    -top    [f i i]  extract posterior distribution of the models, save them as mod96files
                     first argument = lowest log likelyhood value to include (<=0.0, 0.0 means all)
                     second argument = highest model number to include (>=0, 0 means all)
                     third argument = include only one model over "step" (>=1)
                     default {default_extract_llkmin}, {default_extract_limit}, {default_extract_step}
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
--test               testing option
'''.format(
    version=version,
    default_Nworkers=mapkwargs['Nworkers'],
    default_Taskset=mapkwargs['Taskset'],
    default_cmap=default_cmap,
    default_rootnames=default_rootnames,
    default_top_llkmin=default_top_llkmin,
    default_top_limit=default_top_limit,
    default_top_step=default_top_step,
    default_pdf_llkmin=default_pdf_llkmin,
    default_pdf_limit=default_pdf_limit,
    default_pdf_step=default_pdf_step,
    default_extract_llkmin=default_extract_llkmin,
    default_extract_limit=default_extract_limit,
    default_extract_step=default_extract_step,
    default_nchain=default_nchain,
    default_nkeep=default_nkeep,
    default_topstep=default_top_step,
    default_parameterization_list=default_parameterization_list,
    default_parameterization=default_parameterization
    )
# -------------------------------------
example="""
# -------------
# Example usage of HerrMet V{version}
# -------------

# 1/ Data
# get the target dispersion curves, resample it between 0.2-1.5 Hz 
# with 15 samples spaced logarithmically in period domain
# adjust uncertainties to 0.1 in logaritmic domain, 
# overwrite target if exists (_HerrMet.target) 
# and display it
HerrMet --target /path/to/my/data/file.surf96 \\
            -resamp 0.2 1.5 15 plog \\
            -lunc 0.1 \\
            -ot \\
            --display

# >> you may edit one or more target files (e.g. _HerrMet_t???/_HerrMet.target) 
#    and remove points that 
#    do not need to be inverted, check with  
HerrMet --display _HerrMet_t???

# 2/ Parameterization
# build parameter file from existing depthmodel,
# use 4 layers, use parametrization mZVSPRRH, 
# require vp, vs and density to be growing
# overwrite paramfile if exists (_HerrMet_*/_HerrMet.param) and display
HerrMet --param 4 3. \\
            -basedon /path/to/my/depthmodel.mod96 \\
            -t  mZVSPRRH \\
            -growing \\
            -op \\
            --display .

# >> now edit _HerrMet.param and customize it, check with 
HerrMet --display . 

# when ok, send the parameterization file to the rootnames
HerrMet --send

# 3/ Inversion
# run inversion with 4 chains, keep 200 models each, 
# run on 4 virtual threads
HerrMet -w 4 \\
        --run  \\
            -nchain 4 \\
            -nkeep 200

# 4/ Results
# display the best 10 models,
# recompute the disp curves with higher resolution
# compute median and percentiles over the whole population of models
# save as png file (non display backend) 
HerrMet --display \\
            -top 0.0 10 1 \\
            -overdisp \\
            -pdf \\
            -png
""".format(version=version)


# -------------------------------------
def display_function(rootname, argv, verbose):
    targetfile = "%s/_HerrMet.target" % rootname
    paramfile = "%s/_HerrMet.param" % rootname
    runfile = '%s/_HerrMet.run' % rootname
    pngfile = '%s/_HerrMet.png' % rootname

    # ------ Initiate the displayer using the target data if exists
    if os.path.exists(targetfile):
        rd = DepthDispDisplay(targetfile=targetfile)
        d = makedatacoder(targetfile, which=Datacoder_log)  # datacoder based on observations
        dobs, _ = d.target()
    else:
        print "no target file found in %s" % rootname
        rd = DepthDispDisplay()

    # ------ Display run results if exist
    if os.path.exists(runfile) and ("top" in argv.keys() or "pdf" in argv.keys()):

        with RunFile(runfile, verbose=verbose) as rundb:

            # --- display best models
            if "top" in argv.keys():

                assert argv["top"] is None or len(argv["top"]) == 3  # unexpected argument number
                if argv["top"] is None:
                    top_llkmin, top_limit, top_step = default_top_llkmin, default_top_limit, default_top_step
                elif len(argv['top']) == 3:
                    top_llkmin, top_limit, top_step = argv['top']

                print "top : llkmin %f, limit %d, step %d" % (top_llkmin, top_limit, top_step),
                chainids, weights, llks, ms, ds = rundb.like_read_run_1(llkmin=top_llkmin, limit=top_limit,
                                                                        step=top_step)
                vmin, vmax = llks.min(), llks.max()
                colors = values2colors(llks, vmin=vmin, vmax=vmax, cmap=cmap)

                if "overdisp" in argv.keys():
                    """note : recomputing dispersion with another frequency array might
                              result in a completely different dispersion curve in case
                              of root search failure """
                    waves, types, modes, freqs, _ = ds[0]
                    overwaves, overtypes, overmodes, _, _ = zip(
                        *list(groupbywtm(waves, types, modes, freqs, np.arange(len(freqs)), None, True)))
                    overfreqs = [freqspace(0.6 * min(freqs), 1.4 * max(freqs), 100, "plog") for _ in
                                 xrange(len(overwaves))]
                    overwaves, overtypes, overmodes, overfreqs = igroupbywtm(overwaves, overtypes, overmodes, overfreqs)
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

                    cb = makecolorbar(vmin=vmin, vmax=vmax, cmap=cmap)
                    pos = rd.axdisp[-1].get_position()
                    cax = rd.fig.add_axes((pos.x0, 0.12, pos.width, 0.01))
                    rd.fig.colorbar(cb, cax=cax, label="log likelyhood", orientation="horizontal")
                    cax.set_xticklabels(cax.get_xticklabels(), rotation=90., horizontalalignment="center")

                else:
                    "display the dispersion curves as stored in the database"
                    for i in range(len(llks))[::-1]:
                        rd.plotmodel(color=colors[i], alpha=1.0, linewidth=3, *ms[i])
                        rd.plotdisp(color=colors[i], alpha=1.0, linewidth=3, *ds[i])

                    cb = makecolorbar(vmin=vmin, vmax=vmax, cmap=cmap)
                    pos = rd.axdisp[-1].get_position()
                    cax = rd.fig.add_axes((pos.x0, 0.12, pos.width, 0.01))
                    rd.fig.colorbar(cb, cax=cax, label="log likelyhood", orientation="horizontal")
                    cax.set_xticklabels(cax.get_xticklabels(), rotation=90., horizontalalignment="center")

            # ---- display posterior pdf
            if "pdf" in argv.keys():

                assert argv["pdf"] is None or len(argv["pdf"]) == 3  # unexpected argument number
                if argv["pdf"] is None:
                    pdf_llkmin, pdf_limit, pdf_step = default_pdf_llkmin, default_pdf_limit, default_pdf_step
                elif len(argv['pdf']) == 3:
                    pdf_llkmin, pdf_limit, pdf_step = argv['pdf']

                print "pdf : llkmin %f, limit %d, step %d" % (pdf_llkmin, pdf_limit, pdf_step),
                chainids, weights, llks, ms, ds = rundb.like_read_run_1(llkmin=pdf_llkmin, limit=pdf_limit,
                                                                        step=pdf_step)

                dms, wgts = [], []
                for weight, (ztop, vp, vs, rh) in zip(weights, ms):  # readHerrmetout("_HerrMet.run"):
                    dm = depthmodel_from_arrays(ztop, vp, vs, rh)
                    dms.append(dm)
                    wgts.append(weight)
                for p, (vppc, vspc, rhpc, prpc) in \
                        dmstats1(dms,
                                 percentiles=[0.16, 0.5, 0.84],
                                 Ndepth=100,
                                 Nvalue=100,
                                 weights=wgts, **mapkwargs):
                    try:
                        l = 3 if p == 0.5 else 1
                        vppc.show(rd.axvp, color="b", linewidth=l)
                        vspc.show(rd.axvs, color="b", linewidth=l)
                        rhpc.show(rd.axrh, color="b", linewidth=l)
                        prpc.show(rd.axpr, color="b", linewidth=l)

                        # dmout = depthmodel(vppc, vspc, rhpc) #use extract for saveing
                        # dmout.write96('_HerrMet.p%.2f.mod96' % p, overwrite=True)#use extract for saveing
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
    if "m96" in argv.keys():  # plot user data on top
        for m96 in argv['m96']:
            try:
                dm = depthmodel_from_mod96(m96)
                dm.vp.show(rd.axvp, "m", linewidth=3, label=m96)
                dm.vs.show(rd.axvs, "m", linewidth=3)
                dm.rh.show(rd.axrh, "m", linewidth=3)
                dm.pr().show(rd.axpr, "m", linewidth=3)
            except KeyboardInterrupt:
                raise
            except Exception as e:
                print 'could not read or display %s (reason %s)' % (m96, str(e))
            rd.axvp.legend(loc=3)
    # --------------------
    if os.path.exists(targetfile):
        # plot data on top
        rd.plotdisp(d.waves, d.types, d.modes, d.freqs, d.inv(dobs), dvalues=d.dvalues, alpha=0.8, color="g",
                    linewidth=3)

        if "overdisp" in argv.keys():
            rd.set_vlim((0.5 * d.values.min(), 1.5 * d.values.max()))
            rd.set_plim((0.8 / overfreqs.max(), 1.2 / overfreqs.min()))
        else:
            rd.set_vlim((0.8 * d.values.min(), 1.2 * d.values.max()))
            rd.set_plim((0.8 / d.freqs.max(), 1.2 / d.freqs.min()))
    rd.tick()
    rd.grid()
    rd.fig.suptitle(rootname.split('_HerrMet_')[-1])
    chftsz(rd.fig, 10)
    if "png" in argv.keys():
        if verbose:
            print "writing %s" % pngfile
        rd.fig.savefig(pngfile)
    else:
        showme()
    plt.close(rd.fig)


# -------------------------------------
def run_function(rootname, mode, Nchain, Nkeep, argv, verbose):
    runfile = "%s/_HerrMet.run" % rootname
    paramfile = "%s/_HerrMet.param" % rootname
    targetfile = "%s/_HerrMet.target" % rootname

    if mode == "append" and not os.path.exists(runfile):
        mode = "restart"
    elif mode == "restart" and os.path.exists(runfile):
        os.remove(runfile)

    # ------
    p, logRHOM = load_paramfile(paramfile)
    # ------
    d = makedatacoder(targetfile, which=Datacoder_log)  # datacoder based on observations
    dobs, CDinv = d.target()
    duncs = CDinv ** -.5
    ND = len(dobs)
    dinfs = d(0.1 * np.ones_like(d.values))
    dsups = d(3.5 * np.ones_like(d.values))
    logRHOD = LogGaussND(dobs, duncs, dinfs, dsups, k=1000., nanbehavior=1)
    # ------
    G = Theory(parameterizer=p, datacoder=d)
    # ---------------------------------
    if mode == "restart":
        with RunFile(runfile, create=True, verbose=verbose) as rundb:
            rundb.drop()
            rundb.reset(p.NLAYER, d.waves, d.types, d.modes, d.freqs)
    elif mode == "append":
        pass

    # ---------------------------------
    def gen():
        for nchain in xrange(Nchain):
            M0 = np.random.rand(len(p.MINF)) * (p.MSUP - p.MINF) + p.MINF
            MSTD = p.MSTD
            yield Job(nchain, M0, MSTD, nkeep=Nkeep)


    def fun(worker, chainid, M0, MSTD, nkeep):
        models, datas, weights, llks = metropolis(M0, MSTD, G, ND, logRHOD, logRHOM,
              nkeep=nkeep,
              normallaw=worker.randn,
              unilaw=worker.rand,
              chainid=chainid,
              HL=10,
              IK0=0.25,
              MPMIN=1.e-6,
              MPMAX=1e6,
              adjustspeed=0.3,
              nofail=True,
              debug=False,
              verbose=verbose)

        I = np.any(~np.isnan(datas), axis=1)
        models, datas, weights, llks = models[I, :], datas[I, :], weights[I], llks[I]

        return chainid, models, datas, weights, llks
    # ---------------------------------
    with MapAsync(fun, gen(), **mapkwargs) as ma, RunFile(runfile, verbose=verbose) as rundb:
        rundb.begintransaction()

        try:
            for jobid, (chainid, models, datas, weights, llks), _, _ in ma:
                rundb.insert(models, datas, weights, llks, p, d)
                rundb.savepoint()
            rundb.commit()
        except:
            rundb.rollback(crash=True)


# -------------------------------------
def extract_function(rootname, extract_llkmin, extract_limit, extract_step, verbose, percentiles=[.16,.5,.84]):

    percentiles = np.array(percentiles, float)
    assert len(np.unique(percentiles)) == len(percentiles)
    assert np.all(0 < percentiles)
    assert np.all(1 > percentiles)
    assert np.all(percentiles[1:] > percentiles[:-1])

    runfile = "%s/_HerrMet.run" % rootname
    with RunFile(runfile, verbose=verbose) as rundb:
        print "extract : llkmin %f, limit %d, step %d" % (extract_llkmin, extract_limit, extract_step),
        chainids, weights, llks, ms, ds = rundb.like_read_run_1(llkmin=extract_llkmin, limit=extract_limit,
                                                                step=extract_step)

    dms, wgts = [], []
    for weight, (ztop, vp, vs, rh) in zip(weights, ms):  # readHerrmetout("_HerrMet.run"):
        dm = depthmodel(depthmodel1D(ztop, vp),
                        depthmodel1D(ztop, vs),
                        depthmodel1D(ztop, rh))
        dms.append(dm)
        wgts.append(weight)
    for p, (vppc, vspc, rhpc, prpc) in dmstats1(dms,
                                                percentiles=percentiles,
                                                Ndepth=100,
                                                Nvalue=100,
                                                weights=wgts, **mapkwargs):
        try:
            dmout = depthmodel(vppc, vspc, rhpc)
            out = '%s/_HerrMet.p%.2f.mod96' % (rootname, p)
            if verbose:
                print "writing %s" % out
            dmout.write96(out)  # , overwrite=True)
        except KeyboardInterrupt:
            raise
        except Exception as e:
            print "Error", str(e)


# -------------------------------------
if __name__ == "__main__":

    # -------------------------------------
    if len(sys.argv) == 1:
        # no arguments
        print help
        sys.exit()

    # -------------------------------------
    # read arguments
    argv = readargv()
    verbose = default_verbose
    cmap = plt.get_cmap(default_cmap)

    # -------------------------------------
    # prevent typos in arguments, keep the autorizedkeys list up to date
    for k in argv.keys():
        if k not in autorizedkeys:
            if k.startswith('_'): continue # private keys
            raise Exception('keyword %s is not recognized' % k)

    # -------------------------------------
    if "w" in argv.keys():
        mapkwargs["Nworkers"] = int(argv['w'][0])

    # -------------------------------------
    if "taskset" in argv.keys():
        mapkwargs["Taskset"] = argv['taskset'][0]
        os.system('taskset -pc %s %d' % (argv['taskset'][0], os.getpid()))

    # -------------------------------------
    if "lowprio" in argv.keys():
        mapkwargs["LowPriority"] = True

    # -------------------------------------
    if "inline" in argv.keys():
        showme = plt.show

    # -------------------------------------
    if "cmap" in argv.keys():
        try:
            cmap = plt.get_cmap(argv['cmap'][0])
        except ValueError:
            try:
                cmap = eval("cmaps.%s()" % argv['cmap'][0])
            except:
                raise Exception('could not find colormap %s neither in matplotlib neither in tetedenoeud.utils.cmaps' % argv['cmap'][0])

    # -------------------------------------
    if "verbose" in argv.keys():
        if argv['verbose'][0].lower() in ['off', 'false']:
            verbose = False
        elif argv['verbose'][0].lower() in ['on', 'true']:
            verbose = True

    # -------------------------------------
    if "h" in argv.keys() or "help" in argv.keys():
        print help
        sys.exit()

    # -------------------------------------
    if "ex" in argv.keys() or "example" in argv.keys():
        print example
        sys.exit()

    # -------------------------------------
    elif "v" in argv.keys() or "version" in argv.keys():
        print "version : %s" % version
        sys.exit()

    # -------------------------------------
    if "target" in argv.keys():

        # determine root names from target filess
        rootnames = []
        for s96 in argv['target']:
            rootname = "_HerrMet_" + s96.split('/')[-1].split('.s96')[0].split('.surf96')[0]
            rootnames.append(rootname)
        assert len(np.unique(rootnames)) == len(rootnames) # make sure all rootnames are distinct

        # handle already existing files
        if "ot" not in argv.keys():
            for rootname in rootnames:
                if os.path.exists('%s/_HerrMet.target' % rootname):
                    raise Exception('file %s/_HerrMet.target exists already, use -ot to overwrite' % rootname)

        # loop over targets
        for rootname, s96 in zip(rootnames, argv['target']):

            # create directory
            if not os.path.isdir(rootname):
                os.mkdir(rootname)
                os.system('cp %s %s/%s.copy' % (s96, rootname, s96.split('/')[-1]))

            s = surf96reader(s96)
            # -------------------
            if "resamp" in argv.keys():
                news = s.copy()
                news.clear() #clear all entries
                newf = freqspace(freqmin=float(argv["resamp"][0]),
                                 freqmax=float(argv["resamp"][1]),
                                 nfreq=int(argv["resamp"][2]),
                                 scale=argv["resamp"][3])
                for law in s.get_all():
                    law.set_extrapolationmode(1)
                    stdlaw = Claw(freq=law.freq, value=law.dvalue, extrapolationmode=0)

                    newvalues = law(newf)
                    newdvalues = stdlaw(newf)

                    I = ~np.isnan(newvalues)
                    if I.any():
                        N = I.sum()
                        news.data['WAVE'] = np.concatenate((news.data['WAVE'], np.array([law.wave]).repeat(N)))
                        news.data['TYPE'] = np.concatenate((news.data['TYPE'], np.array([law.type]).repeat(N)))
                        news.data['FLAG'] = np.concatenate((news.data['FLAG'], np.array([law.flag]).repeat(N)))
                        news.data['MODE'] = np.concatenate((news.data['MODE'], np.array([law.mode]).repeat(N)))
                        news.data['PERIOD'] = np.concatenate((news.data['PERIOD'], 1. / newf[I]))
                        news.data['VALUE'] = np.concatenate((news.data['VALUE'], newvalues[I]))
                        news.data['DVALUE'] = np.concatenate((news.data['DVALUE'], newdvalues[I]))

                s = news
                # print news
            # -------------------
            if "lunc" in argv.keys():
                # set uncertainties to constant in log domain
                lunc = float(argv["lunc"][0])
                s.data['DVALUE'] = s.data['VALUE'] * lunc
            elif "unc" in argv.keys():
                # set uncertainties to constant in lin domain
                unc = float(argv["unc"][0])
                s.data['DVALUE'] = unc
            # -------------------
            if verbose:
                print "writing %s/_HerrMet.target" % rootname
            s.write96('%s/_HerrMet.target' % rootname)

        # -------------------
        print "please keep only datapoints to invert in */_HerrMet.target"
        print "use option --display to see the target data"
        # sys.exit()

    # -------------------------------------
    if "param" in argv.keys():
        if "op" not in argv.keys():
            assert not os.path.exists('_HerrMet.param')

        nlayer = int(argv["param"][0])
        zbot   = float(argv["param"][1])
        type = argv['t'][0] if "t" in argv.keys() else default_parameterization
        basedon = argv['basedon'][0] if "basedon" in argv.keys() else None
        if "growing" in argv.keys():
            assert "dvs" not in argv.keys() #not compatible with -growing
            assert "dvp" not in argv.keys() #not compatible with -growing
            assert "drh" not in argv.keys() #not compatible with -growing
            assert "dpr" not in argv.keys() # not compatible with -growing
            dvp = 0., 5.
            dvs = 0., 5.
            drh = 0., 5.
            dpr = -5., 0.
        else:
            dvp = minmax(np.asarray(argv['dvp'], float)) if "dvp" in argv.keys() else None
            dvs = minmax(np.asarray(argv['dvs'], float)) if "dvs" in argv.keys() else None
            drh = minmax(np.asarray(argv['drh'], float)) if "drh" in argv.keys() else None
            dpr = minmax(np.asarray(argv['dpr'], float)) if "dpr" in argv.keys() else None

        if not type in default_parameterization_list:
            raise Exception('please pick one type in %s' % str(default_parameterization_list))
        write_default_paramfile(nlayer, zbot, type=type, basedon=basedon, dvp=dvp, dvs=dvs, drh=drh, dpr=dpr)
        print "please customize _HerrMet.param, do not change line orders and metadata"
        print "use option --display to see the depth boundaries"
        # sys.exit()

    # -------------------------------------
    if "send" in argv.keys():
        if not os.path.exists('_HerrMet.param'):
            raise Exception('please use option --param first')

        rootnames = argv['send']
        if rootnames is None:
            rootnames = glob.glob(default_rootnames)
        assert len(rootnames)
        for rootname in rootnames:
            if not os.path.isdir(rootname):
                raise Exception('%s does not exist' % rootname)
            elif not rootname.startswith('_HerrMet_'):
                raise Exception('%s does not starts with _HerrMet_' % rootname)
            elif os.path.exists('%s/_HerrMet.param' % rootname) \
                       and not "op" in argv.keys():
                raise Exception('%s/_HerrMet.param exists already, use -op' % rootname)

        for rootname in rootnames:
            cmd = 'cp ./_HerrMet.param %s/' % rootname
            if verbose:
                print cmd
            os.system(cmd)

    # -------------------------------------
    if "run" in argv.keys():

        mode = "append" if "append" in argv.keys() else "restart"
        rootnames = argv['run']
        if rootnames is None:
            rootnames = glob.glob(default_rootnames)
        assert len(rootnames)

        for rootname in rootnames:
            if not os.path.isdir(rootname):
                raise Exception('%s does not exist' % rootname)
            elif not rootname.startswith('_HerrMet_'):
                raise Exception('%s does not starts with _HerrMet_' % rootname)

        Nchain = int(argv['nchain'][0]) if "nchain" in argv.keys() else default_nchain
        Nkeep = int(argv['nkeep'][0]) if "nkeep" in argv.keys() else default_nkeep

        for rootname in rootnames:
            run_function(rootname, mode=mode, Nchain=Nchain, Nkeep=Nkeep, argv=argv, verbose=verbose)

    # -------------------------------------
    if "extract" in argv.keys():
        rootnames = argv['extract']
        if rootnames is None:
            rootnames = glob.glob(default_rootnames)
        assert len(rootnames)

        for rootname in rootnames:
            if not os.path.isdir(rootname):
                raise Exception('%s does not exist' % rootname)
            elif not rootname.startswith('_HerrMet_'):
                raise Exception('%s does not starts with _HerrMet_' % rootname)

        assert "top" not in argv.keys() or len(argv["top"]) == 3  # unexpected argument number
        if "top" not in argv.keys():
            extract_llkmin, extract_limit, extract_step = default_extract_llkmin, default_extract_limit, default_extract_step
        elif len(argv['top']) == 3:
            extract_llkmin, extract_limit, extract_step = argv['top']

        def gen():
            for rootname in rootnames:
                yield Job(rootname, extract_llkmin, extract_limit, extract_step, verbose)

        with MapAsync(extract_function, gen()) as ma:
            for _ in ma:
                pass

    # -------------------------------------
    if "display" in argv.keys():

        rootnames = argv['display']

        # ----------- special case, just show the parameterization file from --param : ./_HerrMet.param
        if len(rootnames) == 1 and rootnames[0] == '.':
            display_function(".", argv=argv, verbose=verbose)

        # ----------- general case
        else:
            if rootnames is None:
                rootnames = glob.glob(default_rootnames)
            assert len(rootnames)

            for rootname in rootnames:
                if not os.path.isdir(rootname):
                    raise Exception('%s does not exist' % rootname)
                elif not rootname.startswith('_HerrMet_'):
                    raise Exception('%s does not starts with _HerrMet_' % rootname)

            if "png" not in argv.keys():
                # display mode, cannot parallelize
                for rootname in rootnames:
                    display_function(rootname, argv=argv, verbose=verbose)
            else:
                def gen():
                    for rootname in rootnames:
                        yield Job(rootname, argv, verbose=verbose)

                with MapAsync(display_function, gen()) as ma:
                    for _ in ma:
                        pass