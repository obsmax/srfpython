import os, glob
import numpy as np
from srfpython.depthdisp.depthdispdisplay import plt, showme, gcf, gca
from srfpython.HerrMet.files import RunFile

# ------------------------------ defaults
default_rootnames = "_HerrMet_*"

# ------------------------------ autorized_keys
authorized_keys = ["-stats", "-plot", "-delbad", "-delchains"]

# ------------------------------ help messages
short_help = "--manage     summarize run file content, manage run results"

long_help = """\
--manage     s [s..] manage run results for given rootnames, default {default_rootnames}
     -stats          prints detailed stats for each chain of each runfile 
     -plot           display the convergence for every chain and every rootname
     -delbad f       delete bad models, log likelihood below a given threshold, no default
     -delchains i [i...] delete one or more chains using their chainid
      
""".format(default_rootnames=default_rootnames)

# ------------------------------ example usage
example = """\
## MANAGE
# summarize the content of all runfiles,
HerrMet --manage

# print detailed stats and display the convergence 
# of all chains in rootname _HerrMet_001
HerrMet --manage _HerrMet_001 -stats -top

# remove all models whose log likelihood is below -25
# remove chains 8 11 and 14
HerrMet --manage -delbad -25 -delchains 8 11 14
 
"""


# ------------------------------
def manage(argv, verbose, mapkwargs):

    for k in argv.keys():
        if k in ['main', "_keyorder"]:
            continue  # private keys

        if k not in authorized_keys:
            raise Exception('option %s is not recognized' % k)

    rootnames = argv['main']
    if rootnames == []:
        rootnames = glob.glob(default_rootnames)
    assert len(rootnames)
    runfiles = ["%s/_HerrMet.run" % rootname for rootname in rootnames]

    # summarize all files
    for rootname, runfile in zip(rootnames, runfiles):
        with RunFile(runfile, verbose=False) as rundb:
            rundb.summary(head=rootname + " : ")

    # more options
    for rootname, runfile in zip(rootnames, runfiles):

        with RunFile(runfile, verbose=True) as rundb:
            # ------------ print chains stats
            if "-stats" in argv.keys():
                rundb.stats(head=rootname + " : ")

            # ------------ rm
            if "-delbad" in argv.keys():
                rundb.del_bad(llkmin=argv['-delbad'][0])

            if "-delchains" in argv.keys():
                rundb.del_chain(chainid=argv['-delchains'])

            # ------------ plot
            if "-plot" in argv.keys():
                s = rundb.select('''
                select CHAINID, group_concat(NITER), group_concat(LLK) from MODELS
                    group by CHAINID 
                    ''')
                plt.figure()
                for CHAINID, NITER, LLK in s:
                    NITER = np.asarray(NITER.split(','), int)
                    LLK = np.asarray(LLK.split(','), float)
                    gca().plot(NITER, LLK)
                    gca().text(NITER[-1], LLK[-1], CHAINID)
                gca().set_xlabel('# iteration')
                gca().set_ylabel('log likelihood')
                gca().grid(True, linestyle = ":")
                gca().set_title(rootname.split('_HerrMet_')[-1])
                showme()
                plt.close(gcf())
