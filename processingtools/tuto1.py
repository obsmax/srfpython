from multipro7 import *
import time


"""
use keyword arguments
test MapSync
"""

def TargetFunction(x, y = None, *args, **kwargs):
    time.sleep(0.5)
    #do smothing with x, args (tuple) and kwargs (dictionary)
    out = "received : x=%s y=%s args=%s and kwargs=%s" % (x, y, args, kwargs)
    return out

def JobGenerator():
    for jobid in xrange(10):
        x =  2 * jobid
        y =  4 * jobid
        z =  8 * jobid
        t = 10 * jobid
        yield Job(x, y, z, t = t) #use classical syntax for functions, they will be passed to TargetFunction

with MapSync(TargetFunction, JobGenerator()) as ma:
    for jobid, out, _, _ in ma: 
        ma.communicate("job %6d %s" % (jobid, out))

