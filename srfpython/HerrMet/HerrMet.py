#!/usr/bin/env python

"""
HerrMet : code for depth inversion of multimodal/multitypes surface waves
This is the main program to be called from command lines with arguments,
this program will call plugins with corresponding arguments
"""

import os, sys, matplotlib
if "-png" in sys.argv[1:]: matplotlib.use('agg')

from srfpython.utils import readargv1
from srfpython.Herrmann.Herrmann import check_herrmann_codes
from srfpython.HerrMet.plugins import target, param, send, run, \
    manage, neldermead, extract, display, default
check_herrmann_codes()

# -------------------------------------
version = "6.0"
default_verbose=1
default_rootnames = "_HerrMet_*"
default_nworkers = None
default_taskset = None

# -------------------------------------
# security, do not accept options that are not expected => prevent typos
authorized_keys = \
    ["-help", "-h",
     "-version", "-v",
     "-example", "-ex",
     "-w",
     "-taskset",
     "-lowprio",
     "-verbose",
     "--target",
     "--param",
     "--send",
     "--run",
     "--manage",
     "--neldermead",
     "--extract",
     "--display"]

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
{target_help}
{param_help}
{send_help}
{run_help}
{manage_help}
{neldermead_help}
{extract_help}
{display_help}
'''.format(
    version=version,
    default_nworkers=default_nworkers,
    default_taskset=default_taskset,
    default_verbose=default_verbose,
    target_help=target.short_help,
    param_help=param.short_help,
    send_help=send.short_help,
    run_help=run.short_help,
    manage_help=manage.short_help,
    neldermead_help=neldermead.short_help,
    extract_help=extract.short_help,
    display_help=display.short_help)


# -------------------------------------
if __name__ == "__main__":

    argv = readargv1()
    # ------------------------------------- NO ARGUMENT, NAIVE CALL
    if argv == {}:
        # no arguments
        print help
        sys.exit()

    # ------------------------------------- READ ARGUMENTS, CHECK
    # prevent typos in arguments, keep the authorized_keys list up to date
    assert not len(argv['main']) #arguments directly after HerrMet

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

    # ------
    if "--param" in argv.keys():
        param.param(argv['--param'])

    # ------
    if "--send" in argv.keys():
        send.send(argv['--send'], verbose)

    # ------
    if "--run" in argv.keys():
        run.run(argv['--run'], verbose, mapkwargs)

    # ------
    if "--manage" in argv.keys():
        manage.manage(argv['--manage'], verbose, mapkwargs)

    # ------
    if "--neldermead" in argv.keys():
        neldermead.neldermead(argv['--neldermead'], verbose, mapkwargs)

    # ------
    if "--extract" in argv.keys():
        extract.extract(argv['--extract'], verbose, mapkwargs)

    # ------
    if "--display" in argv.keys():
        display.display(argv['--display'], verbose, mapkwargs)
