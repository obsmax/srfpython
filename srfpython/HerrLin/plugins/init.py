import os, glob
from srfpython.depthdisp.depthmodels import depthmodel_from_mod96, depthmodel_from_arrays
from srfpython.HerrMet.files import load_paramfile
from srfpython.HerrLin.files import HerrLinFile

# ------------------------------ defaults
default_rootnames = "_HerrMet_*"

# ------------------------------ autorized_keys
authorized_keys = ["-mean", "-m96", "-oi"]

# ------------------------------ help messages
short_help = "--init       define the starting model for inversion"

long_help = """\
--init       s [s..] set up initial files for inversion, requires list of rootnames to include,
                     default, {default_rootnames}
    -mean            initiate the inversion as the mean of VINF and VSUP columns of the rootnames/_HerrMet.param files
    -m96     s       depth model to be used as initial model (1 file for all rootnames)
    -oi              overwrite rootname/_HerrLin.init if exists   
""".format(default_rootnames=default_rootnames)

# ------------------------------ example usage
example = """\
## INIT
# explain what is done 

HerrLin --init 
"""


# ------------------------------
def init(argv, verbose, mapkwargs):

    # ---------------------
    for k in argv.keys():
        if k in ['main', "_keyorder"]:
            continue  # private keys

        if k not in authorized_keys:
            raise Exception('option %s is not recognized' % k)

    # ---------------------
    rootnames = argv['main']
    if rootnames == []:
        rootnames = glob.glob(default_rootnames)
    assert len(rootnames)

    # ---------------------
    for rootname in rootnames:
        hlffile = "%s/_HerrLin.hlf" % rootname
        if os.path.exists(hlffile):
            if "-oi" not in argv.keys():
                raise Exception('%s exists use -oi to overwrite' % hlffile)
            else:
                os.system('trash %s' % hlffile)

    # ---------------------
    if "-mean" in argv.keys():
        for rootname in rootnames:
            paramfile = "%s/_HerrMet.param" % rootname
            hlffile = "%s/_HerrLin.hlf" % rootname

            p, LogRHOM = load_paramfile(paramfile)
            ztop, vp, vs, rh = p.inv(p.MDEFAULT)
            dm = depthmodel_from_arrays(ztop, vp, vs, rh)

            hlf = HerrLinFile(hlffile)
            if verbose:
                print "init > %s" % hlffile
            hlf.set_init(dm.__str__())

    # ---------------------
    elif "-m96" in argv.keys():
        raise NotImplementedError('problem : the model provided by the user must have at least the same number of layers than the parameterization file (--param --send)')




















        # """initiate all inversions from a given m96 file"""
        # assert len(argv['-m96']) == 1
        #
        # m96 = argv['-m96'][0]
        # assert os.path.exists(m96)
        # dm = depthmodel_from_mod96(m96)
        #
        # for rootname in rootnames:
        #     paramfile = "%s/_HerrMet.param" % rootname
        #     A = AsciiFile(paramfile)
        #     if A.metadata['NLAYER'] != len(dm.vp.z):
        #         print ('Warning : %s is parameterized with %d layers, but you are trying to initiate the inversion with %d layers' % (paramfile, A.metadata['NLAYER'], len(dm.vp.z)))
        #
        #         dm_rootname = dm.copy()
        #
        #     # p, logROM = load_paramfile(paramfile)
        #     #
        #     # assert p.NLAYER == len(dm.vp.z)
        #     #
        #     # b = depthmodel_from_mod96(basedon)
        #     # # ztop = np.linspace(0., b.vp.z.max(), nlayer)
        #     # ztop = depthspace(zbot, nlayer)
        #     # b = b.interp(ztop, interpmethod="weightedstairs")
