#!/usr/bin/env python

#main imports
import sys, matplotlib
if "-agg" in sys.argv[1:] or "--agg" in sys.argv[1:]: matplotlib.use('agg')
import os, imp
import numpy as np

#tetedenoeud
from tetedenoeud.multipro.multipro8 import Job, MapAsync, MapSync
from tetedenoeud.utils.display import values2colors, makecolorbar, legendtext, chftsz
from tetedenoeud.utils.asciifile import AsciiFile

#srfpython
from srfpython.Herrmann.Herrmann import groupbywtm, igroupbywtm
from srfpython.inversion.metropolis2 import LogGaussND, metropolis
from srfpython.depthdisp.depthmodels import depthmodel, depthmodel1D, depthmodel_from_mod96, depthspace
from srfpython.depthdisp.dispcurves import Claw, surf96reader, freqspace
from srfpython.depthdisp.depthpdfs import dmstats1
from srfpython.depthdisp.depthdispdisplay import DepthDispDisplay, plt, gcf, gca, showme, pause
from srfpython.utils import readargv, minmax, tostr

#local imports
from srfpython.HerrMet.files import write_default_paramfile, load_paramfile, read_runfile_1, RunFile
from srfpython.HerrMet.datacoders import makedatacoder, Datacoder_log
from srfpython.HerrMet.theory import Theory, overdisp
from srfpython.HerrMet.parameterizers import Parameterizer_mZVSPRRH, Parameterizer_mZVSVPRH, Parameterizer_mZVSPRzRHz, Parameterizer_mZVSPRzRHvp
#


# -------------------------------------
version = "5"
default_mode = "append"
default_nchain = 12
default_nkeep = 100
default_top = 100
default_topstep = 1
default_parameterization_list = ['mZVSPRRH', 'mZVSVPRH', 'mZVSPRzRHvp', 'mZVSPRzRHz']
default_parameterization = default_parameterization_list[0]
mapkwargs = {"Nworkers" : 4, "Taskset" : "0-3"} #keyword arguments for every parallelized process
# -------------------------------------
autorizedkeys = \
    ["w", "taskset", "agg", "lowprio", "inline",
     "help", "h",
     "example", "ex",
     "param", "basedon", "t", "dvp", "dvs", "drh", "growing", "op",
     "target", "resamp", "lunc", "unc", "ot",
     "run", "nchain", "nkeep", "verbose",
     "extract",
     "disp", "best", "overdisp", "range", "png", "m96",
     "test"]
# -------------------------------------
help = '''HerrMet V{version}
# ----------------------------------------------------
-w           i       set the number of virtual workers to use for all parallelized
                     processes, default {default_Nworkers}
-taskset     s       change job affinity for all parallelized processes, 
                     default {default_Taskset}
-agg                 use agg backend (no display) if mentioned
-lowprio             run processes with low priority if mentioned
-inline              replace showme by plt.show (e.g. jupyter)
-verbose off         turn verbose mode off
# ----------------------------------------------------
--help, -h           display this help message, and quit
--example, -ex       display an example of script, and quit
--param      i f     generate a template parameter file to custom, need the number of layers 
                     and bottom depth in km
    -basedon s       build parametrization based on an existing mod96 file, require a filename, 
                     if not specified, I take fixed values to build the parameter file
    -t       s       parameterization type to use ({default_parameterization_list}), 
                     default {default_parameterization}
                     mZVSPRRH = parameterize with depth interface, 
                                VS in each layer, VP/VS in each layer, Density in each layer
                     mZVSVPRH = parameterize with depth interface,   
                                VP in each layer, VP/VS in each layer, Density in each layer
                     mZVSPRzRHvp = parameterize with depth interface, 
                                use fixed relations for VP/VS versus depth and Density versus VP
                     mZVSPRzRHz = parameterize with depth interface, 
                                use fixed relations for VP/VS versus depth and Density versus depth 
    -dvp     f f     add prior constraint on the vp offset between layers, 
                     requires the extremal values, km/s
    -dvs     f f     add prior constraint on the vs offset between layers, idem
    -drh     f f     add prior constraint on the density offset between layers, idem, g/cm3
    -dpr     f f     add prior constraint on the vp/vs offset between layers, idem, no unit
    -growing         shortcut for -dvp 0. 5. -dvs 0. 5. -drh 0. 5. -dpr -5. 0.
    -op              force overwriting _HerrMet.param if exists
--target     s       set the target dispersion curve from a surf96 file (not modified)
                     the data will be reproduced into a target file that can be customized manually
                     (to remove unwanted points, resample dispersion curves...) 
    -resamp  f f i s resample the dispersion curve in the target file, 
                     requires fmin(Hz), fmax(Hz), nfreq, fscale 
                     (flin=linear in freq domain, plin=linear in period or log=logarithmic scaling)
    -lunc    f       set constant uncertainties in log domain (uncertainty = value x lunc)
    -unc     f       set constant uncertainty in linear domain (uncertainty = unc)
    -ot              force overwriting _HerrMet.target if exists
--run        s       start inversion, requires a running mode 
                     append  : add new models to the exiting run file 
                     restart : overwrite the current run file if any
    -nchain  i       number of chains to use, default {default_nchain}
    -nkeep   i       number of models to keep per chain, default {default_nkeep}
    -w       i       see above, controls the max number of chains to run simultaneously
--extract   [i] [i]  extract posterior distribution on the models, save them as mod96files
                     first argument = number of best models to use/display, default {default_top}
                     second argument = step between them, default {default_topstep}
--disp [i] [i]       display param, target, and run outputs if exists
                     first argument = number of best models to use/display, default {default_top}
                     second argument = step between them, default {default_topstep}
    -best/-overdisp  show the best models on the figure, use overdisp instead of best 
                     to recompute dispersion curves with higher resolution
    -range           compute and show the statistics for the selected models, 
                     use --extract for saving
    -png             save figure as pngfile instead of displaying it on screen
    -m96 file(s)     append depth model(s) to the plot from mod96 file(s)    
--test               testing option
'''.format(
    version=version,
    default_Nworkers=mapkwargs['Nworkers'],
    default_Taskset=mapkwargs['Taskset'],
    default_top=default_top,
    default_nchain=default_nchain,
    default_nkeep=default_nkeep,
    default_topstep=default_topstep,
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
            --disp

# >> you may edit _HerrMet.target and remove points that 
#    do not need to be inverted, check with  
HerrMet --disp

# 2/ Parameterization
# build parameter file from existing depthmodel,
# use 7 layers, use parametrization mZVSPRRH, 
# require vp, vs and density to be growing
# overwrite paramfile if exists (_HerrMet.param) and display
HerrMet --param 7 3. \\
            -basedon /path/to/my/depthmodel.mod96 \\
            -t  mZVSPRRH \\
            -growing \\
            -op \\
            --disp

# >> now edit _HerrMet.param and customize it, check with 
HerrMet --disp

# 3/ Inversion
# run inversion with 12 chains, keep 1000 models each, 
# run on 24 virtual threads
HerrMet -w 24 \\
        --run restart \\
            -nchain 12 -nkeep 1000 

# 4/ Results
# display best 1000 models,
# recompute the best disp curves with higher resolution
# compute median and percentiles over these 1000 models
# save as png file, use a non display backend 
HerrMet -agg \\
        --disp 1000 \\
            -overdisp \\
            -range \\
            -png
""".format(version=version)


# -------------------------------------
if __name__ == "__main__":
    if len(sys.argv) == 1:
        print help
        sys.exit()
    argv = readargv()
    verbose = True
    # nworkers = int(argv['w'][0]) if "w" in argv.keys() else default_nworkers
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
    if "verbose" in argv.keys():
        if argv['verbose'][0].lower() in ['off', 'false']:
            verbose = False
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
    if "clean" in argv.keys():
        os.system('rm -f ./_HerrMet.param ./_HerrMet.target ./_HerrMet.run ./_HerrMet.p*.mod96')
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
        print "use option --disp to see the depth boudaries"
        # sys.exit()
    # -------------------------------------
    if "target" in argv.keys():
        if "ot" not in argv.keys():
            assert not os.path.exists('_HerrMet.target')

        s96 = argv["target"][0]
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
            s.write96('_HerrMet.target')
        elif "unc" in argv.keys():
            # set uncertainties to constant in lin domain
            unc = float(argv["unc"][0])
            s.data['DVALUE'] = unc
        # -------------------
        s.write96('_HerrMet.target')
        print "please only datapoints to invert in _HerrMet.target"
        print "use option --disp to see the target data"
        # sys.exit()
    # -------------------------------------
    if "run" in argv.keys():
        mode = argv['run'][0]
        assert mode in ['append', 'restart']

        Nchain = int(argv['nchain'][0]) if "nchain" in argv.keys() else default_nchain
        Nkeep  = int(argv['nkeep'][0]) if "nkeep" in argv.keys() else default_nkeep

        # ------
        p, logRHOM = load_paramfile('_HerrMet.param')
        # ------
        d = makedatacoder("_HerrMet.target", which=Datacoder_log)  # datacoder based on observations
        dobs, CDinv = d.target()
        duncs = CDinv ** -.5
        ND = len(dobs)
        dinfs = d(0.1 * np.ones_like(d.values))
        dsups = d(3.5 * np.ones_like(d.values))
        logRHOD = LogGaussND(dobs, duncs, dinfs, dsups, k=1000., nanbehavior=1)
        # ------
        G = Theory(parameterizer=p, datacoder=d)
        # ------
        # rd = DepthDispDisplay(targetfile="_HerrMet.target")
        # rd.plotmodel(alpha = 1.0, color = "r", linewidth = 3, *p.inv(p.MINF))
        # rd.plotmodel(alpha = 1.0, color="r", linewidth=3, *p.inv(p.MSUP))
        # rd.plotdisp(d.waves, d.types, d.modes, d.freqs, d.inv(dobs), dvalues = d.dvalues, alpha = 1.0, color="k", linewidth=3)
        # showme(False)
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
        with MapAsync(fun, gen(), **mapkwargs) as ma, open('_HerrMet.run', 'w' if mode == "restart" else "a") as fid:
            fid.write('#CHAINID WEIGHT NLAYER LLK ZTOP[1:] VP VS RH WAVES TYPES MODES FREQS VALUES\n')
            for jobid, (chainid, models, datas, weights, llks), _, _ in ma:
                for mdl, dat, wgt, llk in zip(models, datas, weights, llks):
                    ztop, vp, vs, rh = p.inv(mdl)
                    values = d.inv(dat)
                    nlayer = len(ztop)
                    fid.write("%d %d %d %f %s %s %s %s %s %s %s %s %s\n" %
                              (chainid, wgt, nlayer, llk,
                               tostr(ztop[1:], "%.4f"),
                               tostr(vp, "%.3f"),
                               tostr(vs, "%.3f"),
                               tostr(rh, "%.3f"),
                               tostr(d.waves, "%s"),
                               tostr(d.types, "%s"),
                               tostr(d.modes, "%d"),
                               tostr(d.freqs, "%.4f"),
                               tostr(values, "%4f")))

        # sys.exit()

    # -------------------------------------
    if "test" in argv.keys():
        mode = argv['test'][0]
        assert mode in ['append', 'restart']

        Nchain = int(argv['nchain'][0]) if "nchain" in argv.keys() else default_nchain
        Nkeep  = int(argv['nkeep'][0]) if "nkeep" in argv.keys() else default_nkeep

        # ------
        p, logRHOM = load_paramfile('_HerrMet.param')
        # ------
        d = makedatacoder("_HerrMet.target", which=Datacoder_log)  # datacoder based on observations
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
            with RunFile('_HerrMet.run', create=mode == "restart") as rundb:
                rundb.drop()
                rundb.reset(p.NLAYER, d.waves, d.types, d.modes, d.freqs)
        else:
            raise NotImplementedError('')

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
        with MapAsync(fun, gen(), **mapkwargs) as ma, RunFile("_HerrMet.run") as rundb:
            rundb.begintransaction()

            try:
                for jobid, (chainid, models, datas, weights, llks), _, _ in ma:
                    rundb.insert(models, datas, weights, llks, p, d)
                    rundb.savepoint()
                rundb.commit()
            except:
                rundb.rollback(crash=True)


    # -------------------------------------
    if "extract" in  argv.keys():

        top = int(argv['extract'][0]) if argv['extract'] is not None else default_top
        topstep = int(argv['extract'][1]) if argv['extract'] is not None and len(argv['extract']) >= 2 else default_topstep

        chainids, weights, llks, ms, ds = read_runfile_1('_HerrMet.run', top=top, topstep=topstep)
        dms, wgts = [], []
        for weight, (ztop, vp, vs, rh) in zip(weights, ms):  # readHerrmetout("_HerrMet.run"):
            dm = depthmodel(depthmodel1D(ztop, vp),
                            depthmodel1D(ztop, vs),
                            depthmodel1D(ztop, rh))
            dms.append(dm)
            wgts.append(weight)
        for p, (vppc, vspc, rhpc, prpc) in dmstats1(dms,
                percentiles=[0.16, 0.5, 0.84],
                Ndepth=100,
                Nvalue=100,
                weights=wgts, **mapkwargs):
            try:
                dmout = depthmodel(vppc, vspc, rhpc)
                dmout.write96('_HerrMet.p%.2f.mod96' % p)#, overwrite=True)
            except KeyboardInterrupt: raise
            except Exception as e:
                print "Error", str(e)
    # -------------------------------------
    if "disp" in argv.keys():
        assert not ("best" in argv.keys() and "overdisp" in argv.keys()) #options are not compatible

        top = int(argv['disp'][0]) if argv['disp'] is not None else default_top
        topstep = int(argv['disp'][1]) if argv['disp'] is not None and len(argv['disp']) >= 2 else default_topstep

        # ------ Display the target data if exists
        if os.path.exists("_HerrMet.target"):
            rd = DepthDispDisplay(targetfile="_HerrMet.target")
            d = makedatacoder("./_HerrMet.target", which=Datacoder_log)  # datacoder based on observations
            dobs, _ = d.target()
        else:
            rd = DepthDispDisplay()
            print "call option --target to see the target dispersion curves"

        # ------ Display run results if exist
        if ("best" in argv.keys() or "range" in argv.keys() or "overdisp" in argv.keys()) and os.path.exists('_HerrMet.run'):

            chainids, weights, llks, ms, ds = read_runfile_1('_HerrMet.run', top=top, topstep=topstep)
            vmin, vmax, cmap = llks.min(), llks.max(), plt.cm.gray  #plt.cm.jet# plt.cm.gray #
            colors = values2colors(llks, vmin=vmin, vmax=vmax, cmap=cmap)

            if "best" in argv.keys():
                for i in range(len(llks))[::-1]:
                    ztop, vp, vs, rh = ms[i]
                    dm = depthmodel(depthmodel1D(ztop, vp),
                                    depthmodel1D(ztop, vs),
                                    depthmodel1D(ztop, rh))
                    # dm.write96('M%010.0f.mod' % i, overwrite = True)
                    rd.plotmodel(color=colors[i], alpha=1.0, linewidth=3, *ms[i])
                    rd.plotdisp(color=colors[i], alpha=1.0, linewidth=3, *ds[i])

                cb = makecolorbar(vmin=vmin, vmax=vmax, cmap=cmap)
                cax = rd.fig.add_axes((0.68, 0.1, 0.005, 0.3))
                rd.fig.colorbar(cb, cax=cax, label="log likelyhood")

            elif "overdisp" in argv.keys():
                """note : recomputing dispersion with another frequency array might
                          result in a completely different dispersion curve in case
                          of root search failure """
                waves, types, modes, freqs, _ = ds[0]
                overwaves, overtypes, overmodes, _, _ = zip(*list(groupbywtm(waves, types, modes, freqs, np.arange(len(freqs)), None, True)))
                overfreqs = [freqspace(0.6 * min(freqs), 1.4 * max(freqs), 100, "plog") for _ in xrange(len(overwaves))]
                overwaves, overtypes, overmodes, overfreqs = igroupbywtm(overwaves, overtypes, overmodes, overfreqs)
                for clr, (mms, dds) in zip(colors[::-1], overdisp(ms[::-1], overwaves, overtypes, overmodes, overfreqs,
                                                                  verbose=verbose, **mapkwargs)):
                    rd.plotmodel(color=clr, alpha=1.0, linewidth=3, *mms)
                    try:
                        rd.plotdisp(color=clr, alpha=1.0, linewidth=3, *dds)
                    except KeyboardInterrupt: raise
                    except Exception as e:
                        print "Error : could not plot dispersion curve (%s)" % str(e)

                cb = makecolorbar(vmin=vmin, vmax=vmax, cmap=cmap)
                pos = rd.axdisp[-1].get_position()
                cax = rd.fig.add_axes((pos.x0, 0.1, pos.width, 0.01))
                rd.fig.colorbar(cb, cax=cax, label="log likelyhood", orientation="horizontal")
                cax.set_xticklabels(cax.get_xticklabels(), rotation=90., horizontalalignment="center")

            if "range" in argv.keys():
                dms, wgts = [], []
                for weight, (ztop, vp, vs, rh) in zip(weights, ms):#readHerrmetout("_HerrMet.run"):
                    dm = depthmodel(depthmodel1D(ztop, vp),
                                    depthmodel1D(ztop, vs),
                                    depthmodel1D(ztop, rh))
                    dms.append(dm)
                    wgts.append(weight)
                for p, (vppc, vspc, rhpc, prpc) in dmstats1(dms,
                            percentiles=[0.16, 0.5, 0.84],
                            Ndepth=100,
                            Nvalue=100,
                            weights=wgts, **mapkwargs):
                    # for _ in vppc, vspc, rhpc, prpc:
                    #     _.blur(0.1)
                    try:
                        l = 3 if p == 0.5 else 1
                        vppc.show(rd.axvp, color="b", linewidth=l)
                        vspc.show(rd.axvs, color="b", linewidth=l)
                        rhpc.show(rd.axrh, color="b", linewidth=l)
                        prpc.show(rd.axpr, color="b", linewidth=l)

                        #dmout = depthmodel(vppc, vspc, rhpc) #use extract for saveing
                        #dmout.write96('_HerrMet.p%.2f.mod96' % p, overwrite=True)#use extract for saveing
                    except KeyboardInterrupt: raise
                    except Exception as e:
                        print "Error", str(e)
        else:
            print "call option --run to start inversion"

        # ------
        if os.path.exists('_HerrMet.param'):
            p, _ = load_paramfile('_HerrMet.param')
            showvp, showvs, showrh, showpr = True, True, True, True
            if   isinstance(p, Parameterizer_mZVSVPRH): showpr = False
            elif isinstance(p, Parameterizer_mZVSPRRH): showvp = False
            elif isinstance(p, Parameterizer_mZVSPRzRHvp): showvp = showpr = showrh = False
            elif isinstance(p, Parameterizer_mZVSPRzRHz):  showvp = showpr = showrh = False
            else: raise Exception('')

            # inexact
            # rd.plotmodel(alpha=1.0, color="r", marker = "o--", linewidth=3,
            #              showvp=showvp, showvs=showvs, showrh=showrh, showpr=showpr, *p.inv(p.MINF))
            # rd.plotmodel(alpha=1.0, color="r", marker = "o--", linewidth=3,
            #              showvp=showvp, showvs=showvs, showrh=showrh, showpr=showpr, *p.inv(p.MSUP))
            #
            vplow, vphgh, vslow, vshgh, rhlow, rhhgh, prlow, prhgh = p.boundaries()
            vplow.show(rd.axvp, alpha = 1.0, color = "r", marker = "o--", linewidth = 3)
            vphgh.show(rd.axvp, alpha = 1.0, color = "r", marker = "o--", linewidth = 3)
            vslow.show(rd.axvs, alpha = 1.0, color = "r", marker = "o--", linewidth = 3)
            vshgh.show(rd.axvs, alpha = 1.0, color = "r", marker = "o--", linewidth = 3)
            rhlow.show(rd.axrh, alpha = 1.0, color = "r", marker = "o--", linewidth = 3)
            rhhgh.show(rd.axrh, alpha = 1.0, color = "r", marker = "o--", linewidth = 3)
            prlow.show(rd.axpr, alpha = 1.0, color = "r", marker = "o--", linewidth = 3)
            prhgh.show(rd.axpr, alpha = 1.0, color = "r", marker = "o--", linewidth = 3)
            zmax = 1.1 * p.inv(p.MINF)[0][-1]

            if isinstance(p, Parameterizer_mZVSPRzRHvp):
                rd.axpr.plot(
                    p.PRz(np.linspace(0., zmax, 100)),
                    np.linspace(0., zmax, 100), "r--", linewidth = 3)
                legendtext(rd.axpr, p.PRzName, loc=4)
                legendtext(rd.axrh, p.RHvpName, loc=4)
            elif isinstance(p, Parameterizer_mZVSPRzRHz):
                rd.axpr.plot(
                    p.PRz(np.linspace(0., zmax, 100)),
                    np.linspace(0., zmax, 100), "r--", linewidth = 3)
                rd.axrh.plot(
                    p.RHz(np.linspace(0., zmax, 100)),
                    np.linspace(0., zmax, 100), "r--", linewidth = 3)
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
                except KeyboardInterrupt: raise
                except Exception as e:
                    print 'could not read or display %s (reason %s)' % (m96, str(e))
                rd.axvp.legend(loc=3)
        # --------------------
        if os.path.exists("_HerrMet.target"):
            # plot data on top
            rd.plotdisp(d.waves, d.types, d.modes, d.freqs, d.inv(dobs), dvalues=d.dvalues, alpha=0.8, color="g",linewidth=3)

            if "overdisp" in argv.keys():
                rd.set_vlim((0.5 * d.values.min(), 1.5 * d.values.max()))
                rd.set_plim((0.8 / overfreqs.max(), 1.2 / overfreqs.min()))
            else:
                rd.set_vlim((0.8 * d.values.min(), 1.2 * d.values.max()))
                rd.set_plim((0.8 / d.freqs.max(), 1.2 / d.freqs.min()))
        rd.tick()
        rd.grid()
        chftsz(rd.fig, 10)
        if "png" in argv.keys():
            rd.fig.savefig('_HerrMet.png')
        else:
            showme()
