from __future__ import print_function

import os, glob
import numpy as np
from srfpython.standalone.multipro8 import Job, MapAsync
from srfpython.depthdisp.depthmodels import depthmodel, depthmodel_from_arrays
from srfpython.depthdisp.dispcurves import surf96reader_from_arrays
from srfpython.depthdisp.depthpdfs import dmstats1
from srfpython.depthdisp.disppdfs import dispstats
from srfpython.HerrMet.runfile import RunFile
from srfpython.HerrMet.files import DEFAULTROOTNAMES, HERRMETRUNFILE,\
    HERRMETEXTRACTBESTMODELFILE, HERRMETEXTRACTPDFMODELFILE, HERRMETEXTRACTPDFDISPFILE

# ------------------------------ defaults
default_rootnames = DEFAULTROOTNAMES

default_extract_mode = "last"
default_extract_limit = 0
default_extract_llkmin = 0
default_extract_step = 1

default_top_limit = 10
default_top_llkmin = 0.
default_top_step = 1

# ------------------------------ autorized_keys
authorized_keys = ["-pdf", "-top", "-h", "-help"]

# ------------------------------ help messages
short_help = "--extract    compute and write posterior pdf"

long_help = """\
--extract    s [s..] rootnames for which to compute and extract pdf, default {default_rootnames}
    -pdf   [s i f i] extract posterior distribution of the models, save them as mod96files
                     first argument = selection mode, last or best
                     second argument = highest model number to include (>=0, 0 means all)  
                     third argument = lowest log likelihood value to include (<=0.0, 0.0 means all)
                     fourth argument = include only one model over "step" (>=1)
                     default {default_extract_mode} {default_extract_limit} {default_extract_llkmin} {default_extract_step}
    -top     [i f i] extract best models 
                     first argument = highest model number to include (>=0, 0 means all)
                     second argument = lowest log likelihood value to include (<=0.0, 0.0 means all)  
                     third argument = include only one model over "step" (>=1)
                     default {default_top_limit} {default_top_llkmin} {default_top_step}
    -h, -help        display the help message for this plugin 
""".format(default_rootnames=default_rootnames,
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

HerrMet --extract -pdf best 1000 0. 1 -top 10 0. 1
"""


def _find_runfile(rootname):
    runfile = HERRMETRUNFILE.format(rootname=rootname)
    if not os.path.isfile(runfile):
        err = '{} not found'.format(runfile)
        raise IOError(err)
    return runfile


def _extract_pdf(rootname, extract_mode, extract_limit, extract_llkmin, extract_step, verbose, percentiles, mapkwargs):
    """private"""
    ndepth = 100
    percentiles = np.array(percentiles, float)
    assert len(np.unique(percentiles)) == len(percentiles)
    assert np.all(0 < percentiles)
    assert np.all(1 > percentiles)
    assert np.all(percentiles[1:] > percentiles[:-1])

    for p in percentiles:
        extract_disp_file = HERRMETEXTRACTPDFDISPFILE.format(
            rootname=rootname,
            extract_mode=extract_mode,
            extract_limit=extract_limit,
            extract_llkmin=extract_llkmin,
            extract_step=extract_step,
            percentile=p)
        extract_model_file = HERRMETEXTRACTPDFMODELFILE.format(
            rootname=rootname,
            extract_mode=extract_mode,
            extract_limit=extract_limit,
            extract_llkmin=extract_llkmin,
            extract_step=extract_step,
            percentile=p)
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
        print("extract : llkmin %f, limit %f, step %d" % (extract_llkmin, extract_limit, extract_step))
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
            raise Exception('unexpected extract mode {}'.format(extract_mode))

    if not len(ms):
        return

    dms = [depthmodel_from_arrays(ztop, vp, vs, rh) for ztop, vp, vs, rh in ms]
    for p, (vppc, vspc, rhpc, prpc) in \
            dmstats1(dms,
                     percentiles=percentiles,
                     Ndepth=ndepth,
                     Nvalue=100,
                     weights=weights,
                     **mapkwargs):
        try:
            dmout = depthmodel(vppc, vspc, rhpc)
            extract_model_file = HERRMETEXTRACTPDFMODELFILE.format(
                rootname=rootname,
                extract_mode=extract_mode,
                extract_limit=extract_limit,
                extract_llkmin=extract_llkmin,
                extract_step=extract_step,
                percentile=p)
            if verbose:
                print("writing {}" .format(extract_model_file))
            dmout.write96(extract_model_file)  # , overwrite=True)
        except KeyboardInterrupt:
            raise
        except Exception as e:
            print ("Error", str(e))

    for p in percentiles:
        extract_disp_file = HERRMETEXTRACTPDFDISPFILE.format(
            rootname=rootname,
            extract_mode=extract_mode,
            extract_limit=extract_limit,
            extract_llkmin=extract_llkmin,
            extract_step=extract_step,
            percentile=p)
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
            extract_disp_file = HERRMETEXTRACTPDFDISPFILE.format(
                rootname=rootname,
                extract_mode=extract_mode,
                extract_limit=extract_limit,
                extract_llkmin=extract_llkmin,
                extract_step=extract_step,
                percentile=p)
            if verbose:
                print ("writing to {}".format(extract_disp_file))
            with open(extract_disp_file, 'a') as fid:
                fid.write(srout.__str__())
                fid.write('\n')
        except KeyboardInterrupt:
            raise
        except Exception as e:
            print ("Error", str(e))


def _extract_top(rootname, extract_limit, extract_llkmin, extract_step, verbose):
    """private"""

    try:
        runfile = _find_runfile(rootname)
    except IOError:
        return

    with RunFile(runfile, verbose=verbose) as rundb:
        print ("extract : llkmin %f, limit %d, step %d" % (extract_llkmin, extract_limit, extract_step), end=' ')

        modelids, chainids, weights, llks, nlayers, dms, srs = zip(*list(
            rundb.getpack(limit=extract_limit,
                         llkmin=extract_llkmin,
                         step=extract_step,
                         algo="METROPOLIS")))
        if len(modelids) > 100:
            raise ValueError('too many models to extract')

    if not len(dms):
        return

    for rank, (modelid, chainid, weight, llk, nlayer, dm, sr) in \
            enumerate(zip(modelids, chainids, weights, llks, nlayers, dms, srs)):
        try:
            modelfileout = HERRMETEXTRACTBESTMODELFILE.format(
                rootname=rootname,
                rank=rank,
                modelid=modelid,
                chainid=chainid,
                llk=llk)
            if verbose:
                print ("writing %s" % modelfileout)
            dm.write96(modelfileout)  # , overwrite=True)
        except KeyboardInterrupt:
            raise
        except Exception as e:
            print ("Error", str(e))


def extract(argv, verbose, mapkwargs):

    if '-h' in argv.keys() or "-help" in argv.keys():
        print(long_help)
        return

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
            raise IOError('{} does not exist'.format(rootname))

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
                    percentiles=[0.01, 0.05, 0.16, 0.5, 0.84, 0.95, 0.99], #[.16, .5, .84],
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