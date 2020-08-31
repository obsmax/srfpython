from __future__ import print_function

# ------------------------------ defaults
default_option = None

# ------------------------------ autorized_keys
authorized_keys = ["-option", "-h", "-help"]

# ------------------------------ help messages
short_help = "--default    default plugin structure"

long_help = """\
--default    s [s..] blablabla
    -option  f f i s blablabla, default {default_option}
    -h, -help        display the help message for this plugin 
""".format(default_option=default_option)

# ------------------------------ example usage
example = """\
## DEFAULT
# explain what is done 

HerrMet --default 
"""


def check_keys(argv):
    for k in argv.keys():
        if k in ['main', "_keyorder"]:
            continue  # private keys

        if k not in authorized_keys:
            print('option {} is not recognized'.format(k))  # only for the default plugin

            raise Exception('option {} is not recognized'.format(k))  # please use this line in other plugins


def default(argv, verbose, mapkwargs):

    if '-h' in argv.keys() or "-help" in argv.keys():
        print(long_help)
        return

    check_keys(argv)

    for key, val in argv.items():
        print(key)
        print(val)
        print('')


