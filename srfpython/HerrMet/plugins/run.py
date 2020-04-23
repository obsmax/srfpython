import os, glob
import numpy as np
from srfpython.standalone.multipro8 import Job, MapAsync
from srfpython.inversion.metropolis2 import LogGaussND, metropolis
from srfpython.HerrMet.files import load_paramfile, RunFile
from srfpython.HerrMet.datacoders import makedatacoder, Datacoder_log
from srfpython.HerrMet.theory import Theory

# ------------------------------ defaults
default_rootnames = "_HerrMet_*"
default_nchain = 12
default_nkeep = 100
default_mode = "skip"
# ------------------------------ autorized_keys
authorized_keys = ["-mode", "-nchain", "-nkeep"]

# ------------------------------ help messages
short_help = "--run        invert dispersion data using the Markov Chain Monte Carlo method"

long_help = """\
--run        s [s..] start inversion for the required rootnames, default {default_rootnames}
    -mode    s       set the running mode, default {default_mode}
                     restart  : overwrite the current run file(s) if any   
                     append   : add new models to the exiting run file(s) 
                     skip     : ignore rootnames with existsing run file(s)               
    -nchain  i       number of chains to use, default {default_nchain}
    -nkeep   i       number of models to keep per chain, default {default_nkeep}
    [use -w option before --run to control the maximum number of chains to run simultaneously]   
    """.format(default_rootnames=default_rootnames,
           default_mode=default_mode,
           default_nchain=default_nchain,
           default_nkeep=default_nkeep)

# ------------------------------ example usage
example = """\
## RUN
# run inversion with 4 chains, keep 200 models each, 
# run on 4 virtual threads

HerrMet -w 4 \\
        --run  \\
            -nchain 4 \\
            -nkeep 200
 
"""


# ------------------------------
def run(argv, verbose, mapkwargs):

    for k in argv.keys():
        if k in ['main', "_keyorder"]:
            continue  # private keys

        if k not in authorized_keys:
            raise Exception('option %s is not recognized' % k)

    rootnames = argv['main']
    if rootnames == []:
        rootnames = glob.glob(default_rootnames)
    assert len(rootnames)

    runmode = argv['-mode'][0] if "-mode" in argv.keys() else default_mode
    assert runmode in ['restart', 'append', 'skip']
    Nchain = int(argv['-nchain'][0]) if "-nchain" in argv.keys() else default_nchain
    Nkeep = int(argv['-nkeep'][0]) if "-nkeep" in argv.keys() else default_nkeep

    # ------------------------
    def gen(rootnames, runmode):

        for rootname in rootnames:
            targetfile = "%s/_HerrMet.target" % rootname
            paramfile = "%s/_HerrMet.param" % rootname
            runfile = "%s/_HerrMet.run" % rootname

            if runmode == "append" and not os.path.exists(runfile):
                runmode = "restart"
            elif runmode == "restart" and os.path.exists(runfile):
                os.remove(runfile)
            elif runmode == "skip" and os.path.exists(runfile):
                print "skip %s" % rootname
                continue

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
            if runmode == "restart" or runmode == "skip":
                with RunFile(runfile, create=True, verbose=verbose) as rundb:
                    rundb.drop()
                    rundb.reset(p.NLAYER, d.waves, d.types, d.modes, d.freqs)
            elif runmode == "append":
                pass
            else:
                raise Exception('unexpected runmode %s' % runmode)

            # ---------------------------------
            for chainid in range(Nchain):
                M0 = np.random.rand(len(p.MINF)) * (p.MSUP - p.MINF) + p.MINF
                MSTD = p.MSTD
                yield Job(runfile=runfile,
                          rootname=rootname,
                          chainid=chainid,
                          M0=M0,
                          MSTD=MSTD,
                          G=G,
                          ND=ND,
                          logRHOD=logRHOD,
                          logRHOM=logRHOM,
                          p=p, d=d,
                          nkeep=Nkeep,
                          verbose=verbose)

    # ---------------------------------
    def fun(worker, rootname, runfile,
            chainid, M0, MSTD, G, ND, logRHOD, logRHOM, p, d,
            nkeep, verbose):

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
                                                  verbose=verbose,
                                                  head="%10s " % rootname.split('_HerrMet_')[-1])

        I = np.any(~np.isnan(datas), axis=1)
        models, datas, weights, llks = models[I, :], datas[I, :], weights[I], llks[I]

        return runfile, models, datas, weights, llks, p, d

    # ---------------------------------
    with MapAsync(fun, gen(rootnames, runmode), **mapkwargs) as ma:
        for jobid, answer, _, _ in ma:
            runfile, models, datas, weights, llks, p, d = answer
            if verbose:
                print '=> write to %s' % runfile
            with RunFile(runfile, verbose=False) as rundb:
                rundb.begintransaction()
                try:
                    rundb.insert("METROPOLIS", models, datas, weights, llks, p, d)
                    rundb.commit()
                except:
                    rundb.rollback(crash=True)