import os
import numpy as np
from srfpython.utils import minmax
from srfpython.HerrMet.files import write_default_paramfile

# ------------------------------ defaults
default_parameterization_list = ['mZVSPRRH', 'mZVSVPRH', 'mZVSPRzRHvp', 'mZVSPRzRHz']
default_parameterization = default_parameterization_list[0]

# ------------------------------ autorized_keys
authorized_keys = ["param", "basedon", "t", "dvp", "dvs", "drh", "growing", "op"]

# ------------------------------ help messages
short_help = "--param      create a template parameterization file"

long_help = """\
--param      i f     generate a template parameter file to custom in .
                     need the number of layers and bottom depth in km
    -basedon s       build parametrization based on an existing mod96 file, require a filename, 
                     if not specified, I take fixed values to build the parameter file
    -t       s       parameterization type to use ({default_parameterization_list}), 
                     default {default_parameterization}
          mZVSPRRH = parameterize with 
                     - depth interface (mZ1 = first interface, mZ2, ...), 
                     - VS in each layer (VS0 = first layer, ...), 
                     - VP/VS in each layer (PR0, PR1, ...), 
                     - Density in each layer (RH0, RH1, ...)
          mZVSVPRH = parameterize with depth interface, VP in each layer, 
                     - depth interface (mZ1 = first interface, mZ2, ...), 
                     - VS in each layer (VS0 = first layer, ...), 
                     - VP in each layer (VP0, VP1, ...), 
                     - Density in each layer (RH0, RH1, ...)
       mZVSPRzRHvp = parameterize with depth interface, use fixed 
                     - depth interface (mZ1 = first interface, mZ2, ...), 
                     - VS in each layer (VS0 = first layer, ...), 
                     - use a fixed relation between VP/VS = f(z)
                     - use a fixed relation between RH = f(VP)
        mZVSPRzRHz = parameterize with depth interface, use fixed 
                     - depth interface (mZ1 = first interface, mZ2, ...), 
                     - VS in each layer (VS0 = first layer, ...), 
                     - use a fixed relation between VP/VS = f(z)
                     - use a fixed relation between RH = f(z) 
    -dvp     f f     add prior constraint on the vp offset between layers, 
                     requires the extremal values, km/s
    -dvs     f f     add prior constraint on the vs offset between layers, idem
    -drh     f f     add prior constraint on the density offset between layers, idem, g/cm3
    -dpr     f f     add prior constraint on the vp/vs offset between layers, idem, no unit
    -growing         shortcut for -dvp 0. 5. -dvs 0. 5. -drh 0. 5. -dpr -5. 0.
    -op              force overwriting ./_HerrMet.param if exists
    """.format(default_parameterization_list=default_parameterization_list,
           default_parameterization=default_parameterization)

# ------------------------------ example usage
example = """\
## PARAM
# build parameter file from existing depthmodel,
# use 4 layers, use parametrization mZVSPRRH, 
# require vp, vs and density to be growing
# overwrite paramfile if exists (_HerrMet_*/_HerrMet.param) and display

HerrMet --param 4 3. \\
            -basedon /path/to/my/depthmodel.mod96 \\
            -t  mZVSPRRH \\
            -growing \\
            -op \\
            --display .

# >> now edit _HerrMet.param and customize it, check with 
HerrMet --display . 

# when ok, send the parameterization file to the rootnames
HerrMet --send
"""


# ------------------------------
def param(argv):
    if "op" not in argv.keys():
        assert not os.path.exists('_HerrMet.param')

    nlayer = int(argv["param"][0])
    zbot = float(argv["param"][1])
    type = argv['t'][0] if "t" in argv.keys() else default_parameterization
    basedon = argv['basedon'][0] if "basedon" in argv.keys() else None
    if "growing" in argv.keys():
        assert "dvs" not in argv.keys()  # not compatible with -growing
        assert "dvp" not in argv.keys()  # not compatible with -growing
        assert "drh" not in argv.keys()  # not compatible with -growing
        assert "dpr" not in argv.keys()  # not compatible with -growing
        dvp = 0., 5.
        dvs = 0., 5.
        drh = 0., 5.
        dpr = -5., 0.
    else:
        dvp = minmax(np.asarray(argv['dvp'], float)) if "dvp" in argv.keys() else None
        dvs = minmax(np.asarray(argv['dvs'], float)) if "dvs" in argv.keys() else None
        drh = minmax(np.asarray(argv['drh'], float)) if "drh" in argv.keys() else None
        dpr = minmax(np.asarray(argv['dpr'], float)) if "dpr" in argv.keys() else None

    if not type in default_parameterization_list:
        raise Exception('please pick one type in %s' % str(default_parameterization_list))
    write_default_paramfile(nlayer, zbot, type=type, basedon=basedon, dvp=dvp, dvs=dvs, drh=drh, dpr=dpr)
    print "please customize _HerrMet.param, do not change line orders and metadata"
    print "use option --display to see the depth boundaries"
    # sys.exit()