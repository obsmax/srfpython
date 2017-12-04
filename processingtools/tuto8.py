from multipro7 import *


"""
desactivate paralellization for debogging

simply alias FakeMapAsync as MapAsync in the import section, and comment line when ok
FakeMapAsync has same arguments, attributes and methods than MapSync (and MapAsync) 
but all computations are serial

check CPU activity with htop
"""
from multipro7 import FakeMapAsync as MapAsync #comment me to re-activate parallelization


def TargetFunction(worker):
    t = worker.rand() *4
    start = time.time()
    i = 0
    while time.time() - start < t:
        i + 1 #work for t seconds
    return t
 
def JobGenerator():
    for jobid in xrange(12):
        yield Job()

with MapSync(TargetFunction, JobGenerator()) as ma:
    for n, (jobid, answer, _, _) in enumerate(ma): #
        ma.communicate("job %d returned %s" % (n, answer))

