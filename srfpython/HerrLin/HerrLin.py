#!/usr/bin/env python


import os, sys, matplotlib
if "-png" in sys.argv[1:]: matplotlib.use('agg')

from srfpython.utils import readargv1
from srfpython.Herrmann.Herrmann import check_herrmann_codes
from srfpython.HerrLin.plugins import default, target
check_herrmann_codes()

# -------------------------------------
version = "6.0"
default_verbose=1
default_rootnames = "_HerrLin_*"
default_nworkers = None
default_taskset = None

# -------------------------------------
authorized_keys = \
    ["-help", "-h",
     "-version", "-v",
     "-example", "-ex",
     "-w",
     "-taskset",
     "-lowprio",
     "-verbose",
     "--target"]

# -------------------------------------
help = '''HerrMet V{version}
# ------- main options (s=string, i=int, f=float)
-version, -v          version number, quit
-help, -h   [s...]    help, provide plugin names for details, quit
-example, -ex s [s..] example usage, provide plugin names for details, quit
-w           i        number of virtual workers, default {default_nworkers}
-taskset     s        affinity, e.g. "0-3", default {default_taskset}
-lowprio              run processes with low priority if mentioned
-verbose     i        reduce verbosity, 0/1, default {default_verbose}
# ------- plugins, use --help plugin [plugin ...] for details
{default_help}
'''.format(
    version=version,
    default_nworkers=default_nworkers,
    default_taskset=default_taskset,
    default_verbose=default_verbose,
    default_help=default.short_help)


# -------------------------------------
if __name__ == "__main__":

    argv = readargv1()
    # ------------------------------------- NO ARGUMENT, NAIVE CALL
    if argv == {}:
        # no arguments
        print help
        sys.exit()

    # ------------------------------------- READ ARGUMENTS, CHECK
    # prevent typos in arguments, keep the autorized_keys list up to date
    assert not len(argv['main'])  # arguments directly after HerrMet

    for k in argv.keys():
        if k in ['main', "_keyorder"]:
            continue  # private keys

        if k not in authorized_keys:
            raise Exception('option %s is not recognized' % k)

    # ------------------------------------- VERSION
    if "-v" in argv.keys() or "-version" in argv.keys():
        print "version : %s" % version
        sys.exit()

    # ------------------------------------- HELP
    if "-h" in argv.keys() or "-help" in argv.keys():
        key = "-h" if "-h" in argv.keys() else "-help"
        if argv[key] == []:
            #print full synthteized help
            print help
        else:
            #print specific help for some plugins
            for plugin_name in argv[key]:
                try:
                    print eval('%s.long_help' % plugin_name)
                except NameError:
                    print "%s is not a valid plugin (long_help not found)"
                    continue
                except:
                    raise

        sys.exit()

    # ------------------------------------- EXAMPLES USAGE
    if "-ex" in argv.keys() or "-example" in argv.keys():
        key = "-ex" if "-ex" in argv.keys() else "-example"
        if argv[key] == []:
            raise Exception('please the name of a plugin (e.g. -ex target param send)')
        else:
            #print specific help for some plugins
            for plugin_name in argv[key]:
                try:
                    print eval('%s.example' % plugin_name)
                except NameError:
                    print "%s is not a valid plugin (example not found)"
                    continue
                except:
                    raise

        sys.exit()

        sys.exit()

    # ------------------------------------- MAIN OPTIONS
    mapkwargs = {"Nworkers" : default_nworkers, "Taskset" : default_taskset} #for MapAsync
    verbose = bool(default_verbose)

    # ------
    if "-w" in argv.keys():
        mapkwargs["Nworkers"] = int(argv['-w'][0])

    # ------
    if "-taskset" in argv.keys():
        mapkwargs["Taskset"] = argv['-taskset'][0]
        os.system('taskset -pc %s %d' % (argv['-taskset'][0], os.getpid()))

    # ------
    if "-lowprio" in argv.keys():
        mapkwargs["LowPriority"] = True

    # ------
    if "-verbose" in argv.keys():
        verbose = bool(argv['-verbose'][0])

    # ------------------------------------- PLUGINS
    if "--target" in argv.keys():
        target.target(argv['--target'], verbose)
    #
    # # ------
    # if "param" in argv.keys():
    #     param.param(argv)
    #
    # # ------
    # if "send" in argv.keys():
    #     send.send(argv, verbose)
    #
    # # ------
    # if "run" in argv.keys():
    #     run.run(argv, verbose, mapkwargs)
    #
    # # ------
    # if "extract" in argv.keys():
    #     extract.extract(argv, verbose, mapkwargs)
    #
    # # ------
    # if "display" in argv.keys():
    #     display.display(argv, verbose, mapkwargs)
