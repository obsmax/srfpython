import os, glob
import numpy as np
from srfpython.standalone.multipro8 import Job, MapAsync
from srfpython.depthdisp.depthmodels import depthmodel, depthmodel_from_arrays
from srfpython.depthdisp.dispcurves import surf96reader_from_arrays
from srfpython.depthdisp.depthpdfs import dmstats1
from srfpython.depthdisp.disppdfs import dispstats
from srfpython.HerrMet.files import RunFile


# ------------------------------ defaults
default_rootnames = "_HerrMet_*"
default_extract_mode = "last"
default_extract_limit = 0
default_extract_llkmin = 0
default_extract_step = 1

default_top_limit = 10
default_top_llkmin = 0.
default_top_step = 1

# ------------------------------ autorized_keys
authorized_keys = ["-pdf", "-top"]

# ------------------------------ help messages
short_help = "--extract    compute and write posterior pdf"

long_help = """\
--extract    s [s..] rootnames for which to compute and extract pdf, default {default_rootnames}
    -pdf   [s i f i] extract posterior distribution of the models, save them as mod96files
                     first argument = selection mode, last or best
                     second argument = highest model number to include (>=0, 0 means all)  
                     third argument = lowest log likelyhood value to include (<=0.0, 0.0 means all)
                     fourth argument = include only one model over "step" (>=1)
                     default {default_extract_mode} {default_extract_limit} {default_extract_llkmin} {default_extract_step}
    -top     [i f i] extract best models 
                     first argument = highest model number to include (>=0, 0 means all)
                     second argument = lowest log likelyhood value to include (<=0.0, 0.0 means all)  
                     third argument = include only one model over "step" (>=1)
                     default {default_top_limit} {default_top_llkmin} {default_top_step}
                     """.format(\
           default_rootnames=default_rootnames,
           default_extract_mode=default_extract_mode,
           default_extract_limit=default_extract_limit,
           default_extract_llkmin=default_extract_llkmin,
           default_extract_step=default_extract_step,
           default_top_limit=default_top_limit,
           default_top_llkmin=default_top_llkmin,
           default_top_step=default_top_step)


# ------------------------------ example usage
example = """\
## EXTRACT
# compute pdf using the best 1000 models 

HerrMet --extract -pdf last 1000 0. 1 -top 10 0. 1
"""

RUNFILE = "{rootname}/_HerrMet.run"
EXTRACTMODELFILE = '{rootname}/_HerrMet.p{percentile:.2f}.mod96'
EXTRACTDISPFILE = '{rootname}/_HerrMet.p{percentile:.2f}.surf96'


def _find_runfile(rootname):
    runfile = RUNFILE.format(rootname=rootname)
    if not os.path.isfile(runfile):
        err = '{} not found'.format(runfile)
        raise IOError(err)
    return runfile


def _extract_pdf(rootname, extract_mode, extract_limit, extract_llkmin, extract_step, verbose, percentiles, mapkwargs):
    """private"""
    percentiles = np.array(percentiles, float)
    assert len(np.unique(percentiles)) == len(percentiles)
    assert np.all(0 < percentiles)
    assert np.all(1 > percentiles)
    assert np.all(percentiles[1:] > percentiles[:-1])

    for p in percentiles:
        extract_disp_file = EXTRACTDISPFILE.format(rootname=rootname, percentile=p)
        extract_model_file = EXTRACTMODELFILE.format(rootname=rootname, percentile=p)
        if os.path.isfile(extract_disp_file) and os.path.isfile(extract_model_file):
            continue
        break
    else:
        # did not break => means all files exist already
        print('found existing extraction files in {}, skip'.format(rootname))
        return

    try:
        runfile = _find_runfile(rootname)
    except IOError:
        return

    with RunFile(runfile, verbose=verbose) as rundb:
        print "extract : llkmin %f, limit %d, step %d" % (extract_llkmin, extract_limit, extract_step),
        if extract_mode == "best":
            chainids, weights, llks, ms, ds = \
                rundb.getzip(limit=extract_limit,
                             llkmin=extract_llkmin,
                             step=extract_step,
                             algo="METROPOLIS")
        elif extract_mode == "last":
            chainids, weights, llks, ms, ds = \
                rundb.getlastszip(limit=extract_limit,
                                  llkmin=extract_llkmin,
                                  step=extract_step,
                                  algo="METROPOLIS")
        else:
            raise Exception('unexpected extract mode %s' % extract_mode)

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
            extract_model_file = EXTRACTMODELFILE.format(rootname=rootname, percentile=p)
            if verbose:
                print "writing %s" % extract_model_file
            dmout.write96(extract_model_file)  # , overwrite=True)
        except KeyboardInterrupt:
            raise
        except Exception as e:
            print "Error", str(e)

    for p in percentiles:
        extract_disp_file = EXTRACTDISPFILE.format(rootname=rootname, percentile=p)
        if os.path.isfile(extract_disp_file):
            os.system('trash {}'.format(extract_disp_file))

    for p, (wpc, tpc, mpc, fpc, vpc) in \
            dispstats(ds,
                      percentiles=percentiles,
                      Ndisp=100,
                      weights=weights,
                      **mapkwargs):
        try:
            srout = surf96reader_from_arrays(wpc, tpc, mpc, fpc, vpc)
            extract_disp_file = EXTRACTDISPFILE.format(rootname=rootname, percentile=p)
            if verbose:
                print "writing to {}".format(extract_disp_file)
            with open(extract_disp_file, 'a') as fid:
                fid.write(srout.__str__())
                fid.write('\n')
        except KeyboardInterrupt:
            raise
        except Exception as e:
            print "Error", str(e)


def _extract_top(rootname, extract_limit, extract_llkmin, extract_step, verbose):
    """private"""

    try:
        runfile = _find_runfile(rootname)
    except IOError:
        return

    with RunFile(runfile, verbose=verbose) as rundb:
        print "extract : llkmin %f, limit %d, step %d" % (extract_llkmin, extract_limit, extract_step),

        modelids, chainids, weights, llks, nlayers, dms, srs = zip(*list(
            rundb.getpack(limit=extract_limit,
                         llkmin=extract_llkmin,
                         step=extract_step,
                         algo="METROPOLIS")))
        if len(modelids) > 100:
            raise ValueError('too many models to extract')

    for rank, (modelid, chainid, weight, llk, nlayer, dm, sr) in \
            enumerate(zip(modelids, chainids, weights, llks, nlayers, dms, srs)):
        try:
            out = '{}/_HerrMet.rank{}.modelid{}.chainid{}.llk{}.mod96'.format(
                rootname, rank, modelid, chainid, llk)
            if verbose:
                print "writing %s" % out
            dm.write96(out)  # , overwrite=True)
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
            raise ValueError('option %s is not recognized' % k)

    rootnames = argv['main']
    if not len(rootnames):
        rootnames = glob.glob(default_rootnames)
    assert len(rootnames)

    for rootname in rootnames:
        if not os.path.isdir(rootname):
            raise ValueError('%s does not exist' % rootname)
        elif not rootname.startswith('_HerrMet_'):
            raise ValueError('%s does not start with _HerrMet_' % rootname)

    if "-pdf" in argv.keys():
        if len(argv['-pdf']) == 0:
            extract_mode = default_extract_mode
            extract_limit = default_extract_limit
            extract_llkmin = default_extract_llkmin
            extract_step = default_extract_step

        elif len(argv['-pdf']) == 4:
            extract_mode, extract_limit, extract_llkmin, extract_step = argv['-pdf']
        else:
            raise ValueError('unexpected number of arguments for option -pdf')

        def gen():
            for rootname in rootnames:
                yield Job(
                    rootname, extract_mode, extract_limit,
                    extract_llkmin, extract_step,
                    verbose,
                    percentiles=[.16, .5, .84],
                    mapkwargs=mapkwargs)

        with MapAsync(_extract_pdf, gen(), **mapkwargs) as ma:
            for _ in ma:
                pass

    if "-top" in argv.keys():
        if len(argv['-top']) == 0:
            top_limit, top_llkmin, top_step = \
                default_top_limit, \
                default_top_llkmin, default_top_step
        elif len(argv['-top']) == 3:
            top_limit, top_llkmin, top_step = argv['-top']
        else:
            raise ValueError('unexpected number of arguments for option -top')

        for rootname in rootnames:
            _extract_top(rootname, top_limit, top_llkmin, top_step, verbose)