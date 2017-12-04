from multipro7 import *
from labex.tools.generic_nobspy import primes 
import time

"""
adjust process afinity (i.e. restrict computation on several physical CPUs)
and process priority
adjust number of (virtual) workers
use generator and worker timers
evaluate paralellization efficiency (ratio between processing time and generation time for each job)
"""

def TargetFunction(t):
    "compute prime numbers until time is up"
    start = time.time()
    pgen = primes(1e6)
    while time.time() - start < t:
        p = pgen.next()
    pgen.close()
    return t, p

def JobGenerator():
    for jobid in xrange(100):
        t = np.random.rand() * 2.#pick a random time
        time.sleep(0.1) #slows the generation down, (e.g. data loading)
        yield Job(t)


#use the Taskset argument to adjust the process afinity : here, we use only the first 4 physical CPUs
#the number of workers is virtual and can (slightly) exceed the actual number of CPUs used (e.g. 2 or 3 times more workers than CPUs)

with MapAsync(TargetFunction, JobGenerator(), Taskset="0-3", Nworkers = 8, LowPriority = True) as ma:
    print ma
    for jobid, answer, (gen_start, gen_end), (job_start, job_end) in ma: #
        #gen_start = time at which the generation of the job started
        #gen_end   = time at which the generation of the job ended
        #job_start = time at which the processing of the job started
        #job_end   = time at which the processing of the job ended
        t, p = answer
        adv = (job_end - job_start) / (gen_end - gen_start) 
        ma.communicate("job %6d returned t=%10f p=%10d, job generation took %.8fs, processing took %.8fs, optimal number of workers %d" % \
            (jobid, t, p, gen_end - gen_start, job_end - job_start, adv))

#note that in this example, the processing time (job_end - job_start) corresponds the job argument t as expected
#the processing time must be significantly larger than the generation one
#the higher the ratio between the two durations (adv), the higher the number of workers needed (here we neglect the time used by each iteration of the mapper)
