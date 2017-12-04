from multipro7 import *
import glob, obspy.core as oc
"""
use the asynchronous behavior of the program for file generation 
"""

def filegen():
    """read files in advance and yield them only when needed
       waits if the queue is full
    """
    def jobgen():
        for f in glob.iglob('/data/maximilien/SDS/20*/*/*/???.D/*.*.*.???.D.20*.???'):
            print ">read %s" % f
            st = oc.read(f, format = "MSEED")
            yield Job(f, st)

    def fun(f, st): return f, st #do nothing but pass read data
    with MapSync(fun, jobgen(), Nworkers = 12) as ma:
        for _, (f, st), _, _ in ma:
            yield f, st


print "create file generator"
fg = filegen()
while True:
    if raw_input('next?\n') in ['q', 'quit', 'bye', 'break']: break
    f, st = fg.next()
    print "<%s "% f
fg.close()
