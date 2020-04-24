import os
import numpy as np
from srfpython.depthdisp.dispcurves import Claw, surf96reader, freqspace

# ------------------------------ defaults

# ------------------------------ autorized_keys
authorized_keys = ["-resamp", "-lunc", "-unc", "-ot"]

# ------------------------------ help messages
short_help = "--target     set the target data, create temporary directories for each location"

long_help = """\
--target     s [s..] A target designate a surf96 file containing a set of 
                     dispersion curves to invert for one specific location
                     for each surf96 file, I create a directory in . for temporary files
                     the target file in this directory is a copy of the input surf96 file 
                     it can be edited before inversion
                     (to remove unwanted points or resample dispersion curves...) 
                     it the target file is AAAA.surf96, the temporary directory will be named 
                     _HerrMet_AAAA, this name is referred to as the "rootname" in the later plugins
    -resamp  f f i s resample the dispersion curve in the target file, 
                     requires fmin(Hz), fmax(Hz), nfreq, fscale 
                     (flin=linear in freq domain, plin=linear in period or log=logarithmic scaling)
    -lunc    f       set constant uncertainties in log domain (uncertainty = value x lunc)
    -unc     f       set constant uncertainty in linear domain (uncertainty = unc)
    -ot              force overwriting _HerrMet.target if exists
    """

# ------------------------------ example usage
example = """\
## TARGET
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
            
"""


# ------------------------------
def target(argv, verbose):

    for k in argv.keys():
        if k in ['main', "_keyorder"]:
            continue  # private keys

        if k not in authorized_keys:
            raise Exception('option %s is not recognized' % k)

    # determine root names from target filess
    rootnames = []
    for s96 in argv['main']:
        rootname = "_HerrMet_" + s96.split('/')[-1].split('.s96')[0].split('.surf96')[0]
        rootnames.append(rootname)
    assert len(np.unique(rootnames)) == len(rootnames)  # make sure all rootnames are distinct

    # handle already existing files
    if "-ot" not in argv.keys():
        for rootname in rootnames:
            if os.path.exists('%s/_HerrMet.target' % rootname):
                raise Exception('file %s/_HerrMet.target exists already, use -ot to overwrite' % rootname)

    # loop over targets
    for rootname, s96 in zip(rootnames, argv['main']):

        # create directory
        if not os.path.isdir(rootname):
            os.mkdir(rootname)
            os.system('cp %s %s/%s.copy' % (s96, rootname, s96.split('/')[-1]))

        s = surf96reader(s96)
        # -------------------
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
        # -------------------
        if "-lunc" in argv.keys():
            # set uncertainties to constant in log domain
            lunc = float(argv["-lunc"][0])
            s.data['DVALUE'] = s.data['VALUE'] * lunc
        elif "-unc" in argv.keys():
            # set uncertainties to constant in lin domain
            unc = float(argv["-unc"][0])
            s.data['DVALUE'] = np.ones(len(s.data['VALUE'])) * unc
        # -------------------
        if verbose:
            print "writing %s/_HerrMet.target" % rootname
        s.write96('%s/_HerrMet.target' % rootname)

    # -------------------
    print "please keep only datapoints to invert in */_HerrMet.target"
    print "use option --display to see the target data"