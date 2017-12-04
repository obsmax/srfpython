from multipro7 import *


"""
MapAsync is only suited to cases in which each job should produce a result
StackAsync is more Process-oriented
e.g. you want to stack 100000 objects
     you want 12 workers to stack parts of these objects independently and finally stack the itermediate resutls

the worker may also retain (store) results for each job until the last job has been generated
but this may staturate the internal memory


other application : histogradation 
each worker increments a partial version of the histogram
the last part of the code must merge all partial histograms

"""
#_______________________
class UserStacker(object):
    "the target function designed for internal stack"
    #syntax 1: 
    v = 0. #the initial value of the stack
    #syntax 2:    
    #def __init__(self, v = 0.): 
        #self.v = v #the current value of the stack
    #_______________________
    def __call__(self, worker, x):
        "the way input data (x) should be stacked with current value (self.v)"
        start = time.time()
        while time.time() - start < 0.05: 0 + 0 #slow down, use CPU ressources
        self.v += x
        return self, worker.name #this output will remain in the workspace dedicated to the worker, except for the last one
    #_______________________
    def __iadd__(self, other):
        "a method to merge stackers once back to the serial section"
        assert isinstance(other, UserStacker)
        self.v += other.v
        return self
#_______________________
def JobGen():
    "job generator works like other mappers"
    for i in xrange(1000):
        yield Job(i)
#_______________________
s0 = UserStacker() #the stacker to be reproduced in each independent workspace (deep copy)
with StackAsync(s0, JobGen(), Verbose = False) as sa:
    S = UserStacker() #initiate the final stack
    for jobids, (s, wname), Tgen, Tpro in sa: #receive partial stacks
        sa.communicate("Stacker %s stacked %6d jobs in %.6fs, partial result = %f" % (wname, len(jobids), Tgen + Tpro, s.v))
        S += s #merge all partial stacks using method __iadd__

print "Final result = %f" % S.v
