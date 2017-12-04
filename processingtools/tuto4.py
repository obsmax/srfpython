from multipro7 import *
from labex.tools.generic_nobspy import primes 
import time

"""
use callabe
for instance, you want to attach a large amount of data to the target function
including those data into each job would be an option : i.e.

def gen():
    for jobid in xrange(1000):
        yield Job(data, jobid)
def fun(data, i):
    return data[i]

but this approach implyies to exchange large amount of data between processes and is not optimized

instead, you can attach a copy of the data to each worker
"""

class TargetCallable(object):
    def __init__(self, data):
        self.data = data #attached to the worker
    def __call__(self, i):
        #the function to be parallelized
        #remember that self.data is only a copy of the original data
        #modifying self.data here will have no effect on the original data array
        return self.data[i]

def JobGenerator():
    for jobid in xrange(10):
        yield Job(jobid) 


#create the target object to be attached to each worjer
data = np.arange(1e6) #same data for all jobs
T = TargetCallable(data) #initiate the target function with the data array

with MapAsync(T, JobGenerator()) as ma:
    for jobid, answer, _, _ in ma: #
        data_i = answer
        ma.communicate("job %6d returned %10d" % \
            (jobid, data_i))


#WARNING : a copy of data is attached to each worker, the memory must be large enough to handle this
