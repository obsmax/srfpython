import os

# rootname designates the directory for temporary files of each node
ROOTKEY = "_HerrMet_"
ROOTNAME = os.path.join(".", ROOTKEY + "{node}")

# the parameter file is first custom localy and then sent to the temporary directories
HERRMETPARAMFILELOCAL = "_HerrMet.param"

# the name of the target file in each temporary directory
HERRMETTARGETFILE = os.path.join("{rootname}", "_HerrMet.target")

# the name of the parameter file in each temporary directory
HERRMETPARAMFILE = os.path.join("{rootname}", HERRMETPARAMFILELOCAL)

# the name of the run file in each temporary directory
HERRMETRUNFILE = os.path.join("{rootname}", "_HerrMet.run")

# the name of the display file in each temporary directory
HERRMETDISPLAYFILE = os.path.join("{rootname}", "_HerrMet{options}.png")

# the name of the display file in each temporary directory
HERRMETSTATSFILE = os.path.join("{rootname}", "_HerrMet.stats.png")

# the name of the extraction files in each temporary directory
HERRMETEXTRACTPDFMODELFILE = os.path.join(
    '{rootname}',
    '_HerrMet.{extract_mode:s}_{extract_limit:d}_{extract_llkmin:d}_{extract_step:d}.p{percentile:.2f}.mod96')

HERRMETEXTRACTPDFDISPFILE = os.path.join(
    '{rootname}',
    '_HerrMet.{extract_mode:s}_{extract_limit:d}_{extract_llkmin:d}_{extract_step:d}.p{percentile:.2f}.surf96')

HERRMETEXTRACTBESTMODELFILE = os.path.join(
    '{rootname}',
    '_HerrMet.rank{rank:d}.modelid{modelid:d}.chainid{chainid:d}.llk{llk:f}.mod96')

HERRMETOPTIMIZEPARAMFILE = '_HerrMet.optimize.param'
HERRMETOPTIMIZEPRIORFILE = '_HerrMet.mprior.npz'
HERRMETOPTIMIZEDOBSFILE = '_HerrMet.dobs.npz'
HERRMETOPTIMIZEINVFILE = '_HerrMet.inv.npz'

# the rootnames to use as default in several plugins
DEFAULTROOTNAMES = os.path.dirname(HERRMETTARGETFILE.format(rootname=ROOTNAME.format(node="*")))


def rootname_to_nodename(rootname):
    return os.path.basename(rootname).split(ROOTKEY)[-1]


def surf96filename_to_nodename(surf96filename):
    return os.path.basename(surf96filename).split('.s96')[0].split('.surf96')[0]


def surf96filename_to_rootname(surf96filename):
    return ROOTNAME.format(node=surf96filename_to_nodename(surf96filename))


def surf96filename_to_herrmettargetfile(surf96filename, rootdir="."):
    return HERRMETTARGETFILE.format(
        rootdir=rootdir,
        rootname=surf96filename_to_rootname(surf96filename))


if __name__ == '__main__':
    print(ROOTNAME)
    print(DEFAULTROOTNAMES)
    print(HERRMETPARAMFILE)
    print(HERRMETTARGETFILE)
    print(HERRMETRUNFILE)
    print(HERRMETDISPLAYFILE)
