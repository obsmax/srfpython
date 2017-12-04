from multipro7 import *
from labex.tools.generic_nobspy import primes 
import time

"""
test verbose mode
- verbose can be either simplistic (simply print every operation performed)
- or interactive using curses (once printing as started, use ctrl+H to activate auto-follow mode)
-       job oriented     : displays the status of each job (queued, processed and done)
-       process oriented : displays what is each worker doing.
"""

def TargetFunction(worker):
    t = worker.rand() *4
    if t > 3.: #use communicate to send messages (avoid print for Process or Job verbose modes)
        worker.communicate('Warning (worker %s) : t = %f' % (worker.name, t))
    start = time.time()
    i = 0
    while time.time() - start < t:
        i + 1 #work for t seconds

    return 
 
def JobGenerator():
    for jobid in xrange(100):
        yield Job()

if False: #case 1 : simple verbose mode
    with MapAsync(TargetFunction, JobGenerator(), Verbose=True) as ma:
        for jobid, answer, _, _ in ma: #
            #you may use print here
            ma.communicate('main : received job %d' % jobid)

elif False: #case 2 : Job oriented verbose mode
    with MapAsync(TargetFunction, JobGenerator(), Verbose="Job") as ma:
        for jobid, answer, _, _ in ma: #
            #do not use print!!!
            ma.communicate('main : received job %d' % jobid)

else: #case 3 : Process oriented verbose mode
    with MapAsync(TargetFunction, JobGenerator(), Verbose="Process") as ma:
        for jobid, answer, _, _ in ma: #
            #do not use print!!!
            ma.communicate('main : received job %d' % jobid)
