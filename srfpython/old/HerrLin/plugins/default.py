
# ------------------------------ defaults
default_option = None

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


# ------------------------------
def default(argv, verbose, mapkwargs):

    for k in argv.keys():
        if k in ['main', "_keyorder"]:
            continue  # private keys

        if k not in authorized_keys:
            raise Exception('option %s is not recognized' % k)