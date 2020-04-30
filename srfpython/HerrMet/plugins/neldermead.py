import os, glob
import numpy as np
from srfpython.standalone.multipro8 import Job, MapAsync
from srfpython.inversion.metropolis2 import LogGaussND
from srfpython.inversion.neldermead2 import neldermead as neldermead_function
from srfpython.HerrMet.runfile import RunFile
from srfpython.HerrMet.paramfile import load_paramfile
from srfpython.HerrMet.datacoders import makedatacoder, Datacoder_log
from srfpython.HerrMet.theory import Theory
from srfpython.depthdisp.depthdispdisplay import * #DepthDispDisplay, *

# ------------------------------ defaults
default_rootnames = "_HerrMet_*"
default_top_llkmin = -1000.
default_top_limit = 1
default_top_step = 1
default_niter = 1000
default_ic = 1.e-3
# ------------------------------ autorized_keys
authorized_keys = ["-top", "-niter", "-ic"]

# ------------------------------ help messages
short_help = "--neldermead NOT READY, optimize best solutions from run using the nelder mead algorithm"

long_help = """\
--neldermead s [s..] start inversion for the required rootnames, default {default_rootnames}
        -top    [f i i]  show the best models on the figure, see --extract for arguments
                         default {default_top_llkmin}, {default_top_limit}, {default_top_step}
        -niter  [i]      max number of iteration, default {default_niter}
        -ic     [f]      interruption condition , default {default_ic}        
    """.format(default_rootnames=default_rootnames,
               default_top_llkmin=default_top_llkmin,
               default_top_limit=default_top_limit,
               default_top_step=default_top_step,
               default_niter=default_niter,
               default_ic=default_ic)

# ------------------------------ example usage
example = """\
## NELDERMEAD
# comments 

HerrMet -w 4 \\
        --neldermead
 
"""

def argunique(X):
    return np.unique(X, return_index = True)[1]

# ------------------------------
def neldermead(argv, verbose, mapkwargs):
    raise Exception('not ready')
    for k in argv.keys():
        if k in ['main', "_keyorder"]:
            continue  # private keys

        if k not in authorized_keys:
            raise Exception('option %s is not recognized' % k)

    rootnames = argv['main']
    if rootnames == []:
        rootnames = glob.glob(default_rootnames)
    assert len(rootnames)
    for rootname in rootnames:
        runfile = "%s/_HerrMet.run" % rootname
        assert os.path.exists(runfile)

    assert argv["-top"] == [] or len(argv["-top"]) == 3  # unexpected argument number
    if argv["-top"] == []:
        top_llkmin, top_limit, top_step = default_top_llkmin, default_top_limit, default_top_step
    elif len(argv['-top']) == 3:
        top_llkmin, top_limit, top_step = argv['-top']
    print "top : llkmin %f, limit %d, step %d" % (top_llkmin, top_limit, top_step)

    # ------------------------
    def gen(rootnames):

        for rootname in rootnames:
            targetfile = "%s/_HerrMet.target" % rootname
            paramfile = "%s/_HerrMet.param" % rootname
            runfile = "%s/_HerrMet.run" % rootname

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
            with RunFile(runfile, verbose=verbose) as rundb:
                best = list(rundb.get(llkmin=top_llkmin, limit=top_limit, step=top_step, algo=None))

            # ---------------------------------
            for modelid, chainid, weight, llk, nlayer, model, dat in best:
                M0 = p(*model)
                DM = 1.0 #p.MSTD

                yield Job(runfile=runfile,
                          rootname=rootname,
                          chainid=chainid,
                          M0=M0,
                          DM=DM,
                          G=G,
                          ND=ND,
                          logRHOD=logRHOD,
                          logRHOM=logRHOM,
                          p=p, d=d,
                          verbose=verbose)

    def fun(runfile, rootname, chainid, M0, DM, G, ND, logRHOD, logRHOM, p, d, verbose):
        models, datas, llks = neldermead_function(M0, DM, G, ND, logRHOD, logRHOM,
                                                  alpha=1.0,
                                                  beta=0.9,
                                                  gamma=1.2,
                                                  maxiter=1000,
                                                  interrupt=1e-12,
                                                  debug=1)

        weights = np.ones_like(llks)
        I = np.any(~np.isnan(datas), axis=1)
        models, datas, weights, llks = models[I, :], datas[I, :], weights[I], llks[I]

        I = argunique(llks)
        models, datas, weights, llks = models[I, :], datas[I, :], weights[I], llks[I]

        return runfile, models, datas, weights, llks, p, d


    with MapAsync(fun, gen(rootnames), **mapkwargs) as ma:
        for jobid, answer, _, _ in ma:
            runfile, models, datas, weights, llks, p, d = answer
            if verbose:
                print '=> write to %s' % runfile
            with RunFile(runfile, verbose=False) as rundb:
                rundb.begintransaction()
                try:
                    rundb.insert("NELDERMEAD", models, datas, weights, llks, p, d)
                    rundb.commit()
                except:
                    rundb.rollback(crash=True)


            rd = DepthDispDisplay(targetfile=runfile.replace('.run', '.target'))
            for model, data, llk in zip(models, datas, llks):
                rd.plotmodel(color="k", alpha=0.2, showvp=True, showvs=True, showrh=True,
                          showpr=True, *p.inv(model))

                rd.plotdisp(d.waves,d.types,d.modes,d.freqs,d.inv(data), dvalues=None, color="k", alpha=0.2)
            showme()
            plt.close('all')


