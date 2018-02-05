#!/usr/bin/env python


import os, sys, matplotlib
if "-png" in sys.argv[1:]: matplotlib.use('agg')

from srfpython.utils import readargv
from srfpython.Herrmann.Herrmann import check_herrmann_codes
from srfpython.HerrMet.plugins import target, param, send, run, extract, display, default
check_herrmann_codes()

# -------------------------------------
version = "6.0"
default_verbose=1
default_rootnames = "_HerrMet_*"
default_nworkers = None
default_taskset = None

# -------------------------------------
autorized_keys = \
    ["help", "h",
     "version", "v",
     "example", "ex",
     "w",
     "taskset",
     "lowprio",
     "verbose"]

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
    extract_help=extract.short_help,
    display_help=display.short_help)


# -------------------------------------
if __name__ == "__main__":

    # ------------------------------------- NO ARGUMENT, NAIVE CALL
    if len(sys.argv) == 1:
        # no arguments
        print help
        sys.exit()

    # ------------------------------------- READ ARGUMENTS, CHECK
    argv = readargv()
    # prevent typos in arguments, keep the autorized_keys list up to date
    for k in argv.keys():
        if k.startswith('_'):
            continue  # private keys

        if k not in autorized_keys and \
            k not in target.authorized_keys and \
            k not in param.authorized_keys and \
            k not in send.authorized_keys and \
            k not in run.authorized_keys and \
            k not in extract.authorized_keys and \
            k not in display.authorized_keys:
                raise Exception('option %s is not recognized' % k)

    # ------------------------------------- VERSION
    if "v" in argv.keys() or "version" in argv.keys():
        print "version : %s" % version
        sys.exit()

    # ------------------------------------- HELP
    if "h" in argv.keys() or "help" in argv.keys():
        key = "h" if "h" in argv.keys() else "help"
        if argv[key] is None:
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
    if "ex" in argv.keys() or "example" in argv.keys():
        key = "ex" if "ex" in argv.keys() else "example"
        if argv[key] is None:
            #print full synthteized help
            print help
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
    if "w" in argv.keys():
        mapkwargs["Nworkers"] = int(argv['w'][0])

    # ------
    if "taskset" in argv.keys():
        mapkwargs["Taskset"] = argv['taskset'][0]
        os.system('taskset -pc %s %d' % (argv['taskset'][0], os.getpid()))

    # ------
    if "lowprio" in argv.keys():
        mapkwargs["LowPriority"] = True

    # ------
    if "verbose" in argv.keys():
        verbose = bool(argv['verbose'][0])

    # ------------------------------------- PLUGINS
    if "target" in argv.keys():
        target.target(argv, verbose)

    # ------
    if "param" in argv.keys():
        param.param(argv)

    # ------
    if "send" in argv.keys():
        send.send(argv, verbose)

    # ------
    if "run" in argv.keys():
        run.run(argv, verbose, mapkwargs)

    # ------
    if "extract" in argv.keys():
        extract.extract(argv, verbose, mapkwargs)

    # ------
    if "display" in argv.keys():
        display.display(argv, verbose, mapkwargs)

#     # -------------------------------------
#     if "taskset" in argv.keys():
#         mapkwargs["Taskset"] = argv['taskset'][0]
#         os.system('taskset -pc %s %d' % (argv['taskset'][0], os.getpid()))
#
#     # -------------------------------------
#     if "lowprio" in argv.keys():
#         mapkwargs["LowPriority"] = True
#
#     # -------------------------------------
#     if "verbose" in argv.keys():
#         if argv['verbose'][0].lower() in ['off', 'false']:
#             verbose = False
#         elif argv['verbose'][0].lower() in ['on', 'true']:
#             verbose = True

#
#     # -------------------------------------
#     if "ex" in argv.keys() or "example" in argv.keys():
#         print example
#         sys.exit()
#

#     # -------------------------------------
#     if "target" in argv.keys():
#
#
#
#     # -------------------------------------
#     if "param" in argv.keys():
#
#
#     # -------------------------------------
#     if "send" in argv.keys():
#
#
#     # -------------------------------------
#     if "run" in argv.keys():
#         """"""
#
#
#     # -------------------------------------
#     if "extract" in argv.keys():
#
#
#     # -------------------------------------
#     if "display" in argv.keys():
#
#
#
#
#
#
#
#
#     # # -------------------------------------
#     # if "run" in argv.keys():
#     #
#     #     mode = "append" if "append" in argv.keys() else "restart"
#     #     rootnames = argv['run']
#     #     if rootnames is None:
#     #         rootnames = glob.glob(default_rootnames)
#     #     assert len(rootnames)
#     #
#     #     for rootname in rootnames:
#     #         if not os.path.isdir(rootname):
#     #             raise Exception('%s does not exist' % rootname)
#     #         elif not rootname.startswith('_HerrMet_'):
#     #             raise Exception('%s does not starts with _HerrMet_' % rootname)
#     #
#     #     Nchain = int(argv['nchain'][0]) if "nchain" in argv.keys() else default_nchain
#     #     Nkeep = int(argv['nkeep'][0]) if "nkeep" in argv.keys() else default_nkeep
#     #
#     #     for rootname in rootnames:
#     #         run_function(rootname, mode=mode, Nchain=Nchain, Nkeep=Nkeep, argv=argv, verbose=verbose)
#
# # # -------------------------------------
# # def run_function(rootname, mode, Nchain, Nkeep, argv, verbose):
# #     runfile = "%s/_HerrMet.run" % rootname
# #     paramfile = "%s/_HerrMet.param" % rootname
# #     targetfile = "%s/_HerrMet.target" % rootname
# #
# #     if mode == "append" and not os.path.exists(runfile):
# #         mode = "restart"
# #     elif mode == "restart" and os.path.exists(runfile):
# #         os.remove(runfile)
# #
# #     # ------
# #     p, logRHOM = load_paramfile(paramfile)
# #     # ------
# #     d = makedatacoder(targetfile, which=Datacoder_log)  # datacoder based on observations
# #     dobs, CDinv = d.target()
# #     duncs = CDinv ** -.5
# #     ND = len(dobs)
# #     dinfs = d(0.1 * np.ones_like(d.values))
# #     dsups = d(3.5 * np.ones_like(d.values))
# #     logRHOD = LogGaussND(dobs, duncs, dinfs, dsups, k=1000., nanbehavior=1)
# #     # ------
# #     G = Theory(parameterizer=p, datacoder=d)
# #     # ---------------------------------
# #     if mode == "restart":
# #         with RunFile(runfile, create=True, verbose=verbose) as rundb:
# #             rundb.drop()
# #             rundb.reset(p.NLAYER, d.waves, d.types, d.modes, d.freqs)
# #     elif mode == "append":
# #         pass
# #
# #     # ---------------------------------
# #     def gen():
# #         for nchain in xrange(Nchain):
# #             M0 = np.random.rand(len(p.MINF)) * (p.MSUP - p.MINF) + p.MINF
# #             MSTD = p.MSTD
# #             yield Job(nchain, M0, MSTD, nkeep=Nkeep)
# #
# #
# #     def fun(worker, chainid, M0, MSTD, nkeep):
# #         models, datas, weights, llks = metropolis(M0, MSTD, G, ND, logRHOD, logRHOM,
# #               nkeep=nkeep,
# #               normallaw=worker.randn,
# #               unilaw=worker.rand,
# #               chainid=chainid,
# #               HL=10,
# #               IK0=0.25,
# #               MPMIN=1.e-6,
# #               MPMAX=1e6,
# #               adjustspeed=0.3,
# #               nofail=True,
# #               debug=False,
# #               verbose=verbose)
# #
# #         I = np.any(~np.isnan(datas), axis=1)
# #         models, datas, weights, llks = models[I, :], datas[I, :], weights[I], llks[I]
# #
# #         return chainid, models, datas, weights, llks
# #     # ---------------------------------
# #     with MapAsync(fun, gen(), **mapkwargs) as ma, RunFile(runfile, verbose=verbose) as rundb:
# #         rundb.begintransaction()
# #
# #         try:
# #             for jobid, (chainid, models, datas, weights, llks), _, _ in ma:
# #                 rundb.insert(models, datas, weights, llks, p, d)
# #                 rundb.savepoint()
# #             rundb.commit()
# #         except:
# #            rundb.rollback(crash=True)