from __future__ import print_function

import os
import numpy as np
from srfpython.depthdisp.dispcurves import Claw, surf96reader, freqspace
from srfpython.HerrMet.files import ROOTNAME, HERRMETTARGETFILE, \
    surf96filename_to_rootname, surf96filename_to_herrmettargetfile, \
    surf96filename_to_nodename

# ------------------------------ defaults

# ------------------------------ autorized_keys
authorized_keys = ["-resamp", "-lunc", "-unc", "-ot"]

# ------------------------------ help messages
short_help = "--target     set the target data, create temporary directories for each location"

long_help = """\
--target     s [s..] A target designate a surf96 file containing a set of 
                     dispersion curves to invert for one specific location.
                     For each surf96 file, I create a directory in . for temporary files
                     the target file in this directory is a copy of the input surf96 file. 
                     It can be edited before inversion to remove unwanted points or resample dispersion curves... 
                     it the target file is {surf96file}, the temporary directory will be named 
                     {rootname}, this name is referred to as the "rootname" in the other plugins
    -resamp  f f i s resample the dispersion curve in the target file, 
                     requires fmin(Hz), fmax(Hz), nfreq, fscale 
                     (flin=linear in freq domain, plin=linear in period or log=logarithmic scaling)
    -lunc    f       set constant uncertainties in log domain (uncertainty = value x lunc)
    -unc     f       set constant uncertainty in linear domain (uncertainty = unc)
    -ot              force overwriting the targetfiles if exist
    """.format(surf96file="AAAA.surf96",
               rootname=surf96filename_to_rootname("AAAA.surf96"))

# ------------------------------ example usage
example = """\
## TARGET
# get the target dispersion curves, resample it between 0.2-1.5 Hz 
# with 15 samples spaced logarithmically in period domain
# adjust uncertainties to 0.1 in logaritmic domain, 
# overwrite target files if exist ({herrmettargetfile}) 
# and display it

HerrMet --target {surf96filename} \\
            -resamp 0.2 1.5 15 plog \\
            -lunc 0.1 \\
            -ot \\
            --display
            
# >> you may edit one or more target files (among {herrmettargetfile}) 
#    and remove points that 
#    do not need to be inverted, check with  

HerrMet --display {rootname}
            
""".format(
    surf96filename="/path/to/my/data/node???.surf96",
    herrmettargetfile=surf96filename_to_herrmettargetfile("/path/to/my/data/node???.surf96"),
    rootname=ROOTNAME.format(node="node???"))

def target(argv, verbose):
    for k in argv.keys():
        if k in ['main', "_keyorder"]:
            continue  # private keys

        if k not in authorized_keys:
            raise Exception('option %s is not recognized' % k)

    # determine root names from target filess
    rootnames = []
    for s96 in argv['main']:
        rootname = ROOTNAME.format(node=surf96filename_to_nodename(s96))
        rootnames.append(rootname)
    assert len(np.unique(rootnames)) == len(rootnames)  # make sure all rootnames are distinct

    # handle already existing files
    if "-ot" not in argv.keys():
        for rootname in rootnames:
            target_file = HERRMETTARGETFILE.format(rootname=rootname)
            if os.path.exists(target_file):
                raise Exception('file {} exists already, use -ot to overwrite' % target_file)

    # loop over targets
    for rootname, s96 in zip(rootnames, argv['main']):
        target_file = HERRMETTARGETFILE.format(rootname=rootname)
        target_dir = os.path.dirname(target_file)

        # create directory
        if not os.path.isdir(target_dir):
            if verbose:
                print("creating " + target_dir)
            os.makedirs(target_dir)
            s96copy = os.path.join(target_dir, os.path.basename(s96)) + ".copy"
            cpcmd = 'cp {} {}'.format(s96, s96copy)
            if verbose:
                print(cpcmd)
            os.system(cpcmd)

        s = surf96reader(s96)

        if "-resamp" in argv.keys():
            news = s.copy()
            news.clear()  # clear all entries
            newf = freqspace(freqmin=float(argv["-resamp"][0]),
                             freqmax=float(argv["-resamp"][1]),
                             nfreq=int(argv["-resamp"][2]),
                             scale=argv["-resamp"][3])
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

        if "-lunc" in argv.keys():
            # set uncertainties to constant in log domain
            lunc = float(argv["-lunc"][0])
            s.data['DVALUE'] = s.data['VALUE'] * lunc

        elif "-unc" in argv.keys():
            # set uncertainties to constant in lin domain
            unc = float(argv["-unc"][0])
            s.data['DVALUE'] = np.ones(len(s.data['VALUE'])) * unc

        if verbose:
            print("writing " + target_file)
        s.write96(target_file)

    # -------------------
    print("please keep only datapoints to invert in " + HERRMETTARGETFILE.format(rootname="*"))
    print("use option --display to see the target data")
