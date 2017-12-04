from multipro7 import *
from labex.tools.generic_nobspy import primes 
import time

"""
handle errors

- rememeber that an error in the job generator will interrupt the whole process
- you can choose wether an error encountered by one of the workers should be ignored or interrupt the whole process
- an error in the main process (e.g. ctrl+c) will also interrupt the whole process (thanks to the with statement)
"""

def TargetFunction(t):
    if t >= 1.:
        raise ValueError('unexpected value for t')

    start = time.time()
    pgen = primes(1e6)
    while time.time() - start < t:
        p = pgen.next()
    pgen.close()
    return t, p

def JobGenerator():
    for jobid in xrange(10):
        t = np.random.rand() * 2.#pick a random time
##un-comment these two lines to test error encoutered inside the generator
#        if jobid == 3:
#            raise Exception('I am the generator and I don t want to generate that job')
        yield Job(t) 





#case 1 : ignore errors, this case is not compatible with MapSync
print "case 1"
with MapAsync(TargetFunction, JobGenerator(), RaiseIfError = False) as ma:
    for jobid, answer, _, _ in ma: #
        t, p = answer
        ma.communicate("job %6d returned %10f %10d" % \
            (jobid, t, p))
print "case 1 done"





#case 2 : interrupt the whole process if one of the workers raises an exception
print "case 2"
# you can see the error message printed by the worker himself and the the interruption message
with MapAsync(TargetFunction, JobGenerator(), RaiseIfError = True) as ma:
    for jobid, answer, _, _ in ma: #
        t, p = answer
        ma.communicate("job %6d returned %10f %10d" % \
            (jobid, t, p))
print "case 2 done"
