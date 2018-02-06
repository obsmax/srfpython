import os, sys
import numpy as np
from srfpython.HerrLin.files import HerrLinFile
from srfpython.depthdisp.dispcurves import surf96reader

# ------------------------------ defaults

# ------------------------------ autorized_keys
authorized_keys = ["option"]

# ------------------------------ help messages
short_help = "--target     short help"

long_help = """\
--target     s [s..] list of target data (files at surf96 format) to be inverted, 
                     will reproduce file contents into .hlf files with same root names in current directory
                     existing hlf files will be overwriten
    -unc     f       overwrite uncertainties and set it to constant (in km/s, same for all points in the s96 files)
    -lunc    f       overwrite uncertainties and set it to constant in log domain (in km/s)    
"""

# ------------------------------ example usage
example = """\
## TARGET
# explain what is done 

HerrMet --target 
"""


# ------------------------------
def target(argv, verbose):

    for k in argv.keys():
        if k in ['main', "_keyorder"]:
            continue  # private keys

        if k not in authorized_keys:
            raise Exception('option %s is not recognized' % k)

    s96s = argv['main']
    hlfs = []
    for s96 in s96s:
        assert s96.endswith(".surf96")
        assert os.path.exists(s96)
        hlf = s96.split('/')[-1].replace('.surf96', '.hlf')
        hlfs.append(hlf)
    assert len(hlfs) == len(np.unique(hlfs))  # hlf files must have distinct names

    for s96, hlf in zip(s96s, hlfs):
        if os.path.exists(hlf):
            os.system('rm -f %s' % hlf)
        print "%s > %s" % (s96, hlf)

        # -----------------
        s = surf96reader(s96)
        if "-unc" in argv.keys():
            s.data['DVALUE'] = float(argv['-unc'][0]) * np.ones_like(s.data['VALUE'])  # set constant uncertainty
        elif "-lunc" in argv.keys():
            s.data['DVALUE'] = float(
                argv['-lunc'][0]) * s.data.VALUE  # set constant uncertainty to the logarithm of the value

        # -----------------
        h = HerrLinFile(hlf)
        h.set_target(filename=s96, data=s.__str__())
