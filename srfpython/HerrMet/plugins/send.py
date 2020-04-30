from __future__ import print_function
from srfpython.HerrMet.files import HERRMETPARAMFILELOCAL, HERRMETPARAMFILE, DEFAULTROOTNAMES
import os, glob

# ------------------------------ defaults
default_rootnames = DEFAULTROOTNAMES

# ------------------------------ autorized_keys
authorized_keys = ["-op"]

# ------------------------------ help messages
short_help = "--send       send the parameterization to the temporary directories"

long_help = """\
--send       s [s..] send the custom parameterization file {herrmetparamfilelocal} to the specified rootnames, 
                     default {default_rootnames}
    -op              force overwriting {herrmetparamfile} if exists
    """.format(
    herrmetparamfilelocal=HERRMETPARAMFILELOCAL,
    herrmetparamfile=HERRMETPARAMFILE.format(rootname='[rootname]'),
    default_rootnames=default_rootnames)

# ------------------------------ example usage
example = ""


# ------------------------------
def send(argv, verbose):
    for k in argv.keys():
        if k in ['main', "_keyorder"]:
            continue  # private keys

        if k not in authorized_keys:
            raise Exception('option {} is not recognized'.format(k))

    if not os.path.exists(HERRMETPARAMFILELOCAL):
        raise Exception('please use option --param first')

    rootnames = argv['main']
    if rootnames == []:
        rootnames = glob.glob(default_rootnames)
    assert len(rootnames)

    for rootname in rootnames:
        paramfile = HERRMETPARAMFILE.format(rootname=rootname)
        paramdir = os.path.dirname(paramfile)

        if not os.path.isdir(paramdir):
            raise IOError(paramdir + ' is not a directory')

        elif os.path.isfile(paramfile):
            if not "-op" in argv.keys():
                raise Exception('{} exists already, use -op'.format(paramfile))

    for rootname in rootnames:
        paramfile = HERRMETPARAMFILE.format(rootname=rootname)

        cmd = 'cp {} {}'.format(HERRMETPARAMFILELOCAL, paramfile)
        if verbose:
            print(cmd)
        os.system(cmd)