import os, glob
from srfpython.HerrMet.files import load_paramfile
from srfpython.HerrMet.datacoders import makedatacoder, Datacoder_log
from srfpython.HerrMet.theory import Theory
from srfpython.inversion.metropolis2 import LogGaussND
from srfpython.inversion.linearized6 import LinInv
import numpy as np

# ------------------------------ defaults
default_option = None
default_damping = 1.0

# ------------------------------ autorized_keys
authorized_keys = ["-option"]

# ------------------------------ help messages
short_help = "--default    default plugin structure"

long_help = """\
--default    s [s..] blablabla
    -option  f f i s blablabla, default {default_option}
""".format(default_option = default_option)

# ------------------------------ example usage
example = """\
## DEFAULT
# explain what is done 

HerrLin --default 
"""


def _run(rootname, argv, verbose):
    targetfile = "%s/_HerrMet.target" % rootname
    paramfile = "%s/_HerrMet.param" % rootname
    hlffile = "%s/_HerrLin.hlf" % rootname

    damping = argv['damp'][0] if 'damp' in argv.keys() else default_damping

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

    LI = LinInv(g=G,
                d=dobs,
                CDinv=CDinv,
                mapr=,
                CMinv=,
                damping = damping,
                locked = None,
                pertu = 0.05,
                eps = 1.e-6,
                method = 1)
# ------------------------------
def run(argv, verbose, mapkwargs):


    # -----------------
    for k in argv.keys():
        if k in ['main', "_keyorder"]:
            continue  # private keys

        if k not in authorized_keys:
            raise Exception('option %s is not recognized' % k)

    # -----------------
