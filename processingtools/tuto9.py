from labex.processingtools.multipro7 import *
import time


"""

the faster the job generator relative to the processing the better, 
thus parallelization may be usefull inside the job generator
which requires this type of imbricated parallelization scheme

the two levels can work on the same physical CPUs if needed

"""


def TargetFunctionLevel1(t, i):
    start = time.time()
    while time.time() - start < t:
        i += 1 #work for t seconds
    return t, i

def JobGenLevel1():
    #prepare deepest level
    def TargetFunctionLevel0(t):
        start = time.time()
        i = 0
        while time.time() - start < t:
            i += 1 #work for t seconds
        return t, i

    def JobGenLevel0():
        for i in xrange(100):
            t = np.random.rand()
            yield Job(t)

    #start job generation for level1
    with MapAsync(TargetFunctionLevel0, JobGenLevel0(), Taskset="0-11") as ma:
        for jobid, (t, i), _, _ in ma:
            ma.communicate("level0 : job %5d returned %f %d" % (jobid, t, i))
            yield Job(t, i = i) #yields job for level 1 here !!!


#run highest parallelization level
with MapAsync(TargetFunctionLevel1, JobGenLevel1(), Taskset="0-11") as ma:
    for jobid, (t, i), _, _ in ma:
            ma.communicate("level1 : job %5d returned %f %d" % (jobid, t, i))

