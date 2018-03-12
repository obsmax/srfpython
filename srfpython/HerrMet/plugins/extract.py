import os, glob
import numpy as np
from srfpython.standalone.multipro8 import Job, MapAsync
from srfpython.depthdisp.depthmodels import depthmodel, depthmodel_from_arrays
from srfpython.depthdisp.dispcurves import surf96reader_from_arrays
from srfpython.depthdisp.depthpdfs import dmstats1
from srfpython.depthdisp.disppdfs import dispstats
from srfpython.HerrMet.files import  RunFile


# ------------------------------ defaults
default_rootnames = "_HerrMet_*"
default_extract_llkmin = 0
default_extract_limit = 0
default_extract_step = 1

# ------------------------------ autorized_keys
authorized_keys = ["-top"]

# ------------------------------ help messages
short_help = "--extract    compute and write posterior pdf"

long_help = """\
--extract    s [s..] rootnames for which to compute and extract pdf, default {default_rootnames}
    -top    [f i i]  extract posterior distribution of the models, save them as mod96files
                     first argument = lowest log likelyhood value to include (<=0.0, 0.0 means all)
                     second argument = highest model number to include (>=0, 0 means all)
                     third argument = include only one model over "step" (>=1)
                     default {default_extract_llkmin}, {default_extract_limit}, {default_extract_step}
                     """.format(default_rootnames=default_rootnames,
           default_extract_llkmin=default_extract_llkmin,
           default_extract_limit=default_extract_limit,
           default_extract_step=default_extract_step)


# ------------------------------ example usage
example = """\
## EXTRACT
# compute pdf using the best 1000 models 

HerrMet --extract -top 0. 1000 1
"""


# -------------------------------------
def _extract_function(rootname, extract_llkmin, extract_limit, extract_step, verbose, percentiles, mapkwargs):
    """private"""
    percentiles = np.array(percentiles, float)
    assert len(np.unique(percentiles)) == len(percentiles)
    assert np.all(0 < percentiles)
    assert np.all(1 > percentiles)
    assert np.all(percentiles[1:] > percentiles[:-1])

    runfile = "%s/_HerrMet.run" % rootname
    with RunFile(runfile, verbose=verbose) as rundb:
        print "extract : llkmin %f, limit %d, step %d" % (extract_llkmin, extract_limit, extract_step),
        chainids, weights, llks, ms, ds = rundb.getzip(llkmin=extract_llkmin, limit=extract_limit,
                                                                step=extract_step)

    dms = [depthmodel_from_arrays(ztop, vp, vs, rh) for ztop, vp, vs, rh in ms]
    for p, (vppc, vspc, rhpc, prpc) in \
            dmstats1(dms,
                     percentiles=percentiles,
                     Ndepth=100,
                     Nvalue=100,
                     weights=weights,
                     **mapkwargs):
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


    # display the disp pdf
    for p in percentiles:
        out = '%s/_HerrMet.p%.2f.surf96' % (rootname, p)
        os.system('trash %s' % out)
    for p, (wpc, tpc, mpc, fpc, vpc) in \
            dispstats(ds,
                      percentiles=percentiles,
                      Ndisp=100,
                      weights=weights,
                      **mapkwargs):
        try:
            srout = surf96reader_from_arrays(wpc, tpc, mpc, fpc, vpc)
            out = '%s/_HerrMet.p%.2f.surf96' % (rootname, p)
            if verbose:
                print "writing to %s" % out
            with open(out, 'a') as fid:
                fid.write(srout.__str__())
                fid.write('\n')
        except KeyboardInterrupt:
            raise
        except Exception as e:
            print "Error", str(e)


# ------------------------------
def extract(argv, verbose, mapkwargs):
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
        if not os.path.isdir(rootname):
            raise Exception('%s does not exist' % rootname)
        elif not rootname.startswith('_HerrMet_'):
            raise Exception('%s does not starts with _HerrMet_' % rootname)

    assert "-top" not in argv.keys() or len(argv["-top"]) == 3  # unexpected argument number
    if "-top" not in argv.keys():
        extract_llkmin, extract_limit, extract_step = default_extract_llkmin, default_extract_limit, default_extract_step
    elif len(argv['-top']) == 3:
        extract_llkmin, extract_limit, extract_step = argv['-top']

    def gen():
        for rootname in rootnames:
            yield Job(rootname, extract_llkmin, extract_limit, extract_step, verbose,
                      percentiles=[.16,.5,.84], mapkwargs=mapkwargs)

    with MapAsync(_extract_function, gen(), **mapkwargs) as ma:
        for _ in ma:
            pass