from multipro7 import *
from labex.tools.generic_nobspy import primes 
import time

"""
use method and attributes of the worker inside the target function
this is usefull to
- determine the name (id) of the worker that is doing the job
- most important use is for safe random operations (each worker has a different seed)
"""

def TargetFunction(worker, *args, **kwargs):
    'keyword "worker" must appear as first argument'
    sleeptime = worker.rand() #using numpy.random.rand instead would be thread unsafe
    time.sleep(sleeptime) 
    #do something with args and kwargs
    return worker.name, sleeptime, args, kwargs
 
def JobGenerator():
    for jobid in xrange(100):
        yield Job(jobid, "abcdefghijkl"[jobid % 12], k = jobid % 2) 


#create the target object to be attached to each worjer
with MapAsync(TargetFunction, JobGenerator(), Nworkers = 5) as ma:
    for jobid, answer, _, _ in ma: #
        wname, sleeptime, args, kwargs = answer
        ma.communicate('job %d processed by worker named "%s" slept for %fs, got arguments %s and keyword arguments %s' % (jobid, wname, sleeptime, args, kwargs))



#note that several jobs are processed by the same worker !
#The same operation can be done with callables, 
#the "worker" argument must appear after the self argument like:
#class TargetCallable(object):
#   def __init__(self, ...): ...
#   def __call__(self, worker, *args, **kwargs): ...

