import os, glob

# ------------------------------ defaults
default_rootnames = "_HerrMet_*"

# ------------------------------ autorized_keys
authorized_keys = ["send", "op"]

# ------------------------------ help messages
short_help = "--send       send the parameterization to the temporary directories"

long_help = """\
--send       s [s..] send the custom parameterization file ./_HerrMet.param to the specified rootnames, 
                     default {default_rootnames}
    -op              force overwriting ./rootname/_HerrMet.param if exists
    """.format(default_rootnames=default_rootnames)

# ------------------------------ example usage
example = ""


# ------------------------------
def send(argv, verbose):
    if not os.path.exists('_HerrMet.param'):
        raise Exception('please use option --param first')

    rootnames = argv['send']
    if rootnames is None:
        rootnames = glob.glob(default_rootnames)
    assert len(rootnames)
    for rootname in rootnames:
        if not os.path.isdir(rootname):
            raise Exception('%s does not exist' % rootname)
        elif not rootname.startswith('_HerrMet_'):
            raise Exception('%s does not starts with _HerrMet_' % rootname)
        elif os.path.exists('%s/_HerrMet.param' % rootname) \
                and not "op" in argv.keys():
            raise Exception('%s/_HerrMet.param exists already, use -op' % rootname)

    for rootname in rootnames:
        cmd = 'cp ./_HerrMet.param %s/' % rootname
        if verbose:
            print cmd
        os.system(cmd)