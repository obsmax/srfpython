from multipro7 import *
from labex.tools.stdout import waitbar
import time

"""
use waiting bar
"""

def TargetFunction(worker):
    t = worker.rand() *4
    start = time.time()
    i = 0
    while time.time() - start < t:
        i + 1 #work for t seconds
    return 
 
def JobGenerator():
    for jobid in xrange(100):
        yield Job()

with MapAsync(TargetFunction, JobGenerator(), Verbose=False) as ma:
    w = waitbar('progress') #do not print or communicate...
    for n, (jobid, answer, _, _) in enumerate(ma): #
        w.refresh(n / 100.)
    w.close()
