from stdout import InteractiveStdOut
from multiprocessing import queues, Process, Queue, cpu_count
from math import erf as merf
import numpy as np
import copy
import inspect
import os
import random
import sys
import time
import traceback
import types


"""
ML 06-2016 (version 5)
new feature : measures the generation and processing times for each job
can be used to plot the job processing's history or to estimate the optimum worker number to use
Job1 has been renamed has Job!!!

ML 08-2016 (version 5.1)
new feature : the Verbose mode has been improved, messages are redirected to a specific queue
depending on which printer is used, the messages are ignored, simply printed as they come (time ordered, this is the safest verbose mode), 
printed in process-oriented display mode or printed in a job-oriented display mode
one more process is now dedicated to the printing (see InteractiveStdOut)
(there was already one for the input queue and one per worker)
messages in the queue are always with the form : (sender, time, message, jobid), the printer defines how to interpret them

ML 31-08-2016 (version 6)
new feature : 
1) It is now possible to use a callable object as a target (instead of a function)
   This must be used to pack common variables into the target only once; instead of passing the same variables to the target through all jobs
   Remember not to try to modify the target attributes inside the __call__ method, because __call__ is run in a separate workspace
2) You can now use one target per worker to make the worker do different things
   (e.g., each worker may be a markov chain with specific parameters)

ML 12-09-2017 (version 7)
new feature : two more instances : StackAsync (and private WorkerStacker)
StackAsync works like MapAsync except that the worker will run alone without puting results into the output queue until it gots the ending signals
Use it for stacking processes
each worker will internally stack Jobs (no matter which worker stacks what with what)
you get one result per worker with 
jobids = the list of jobs stacked by the worker
answer = the result of the stacking process
Tgen   = the total time spent in jog generation (s)
Tpro   = the total time spent in stacking (s)
=> not tested (RaiseIfError? Taskset? Verbose?, Efficiency?)

ML 28-11-2017
new feature : add the possibility to run processes with low priority


ML 06-01-2018
move multipro to a new independent package (-), rename it as multipro8

ML 17-02-1028
create a local copy of this file for srfpython 
"""


def ut(t):
    s = time.ctime(t)
    return s.split()[3]


def erf(x):
    return np.asarray([merf(w) for w in x], float)


# __________________________________
# Few exception aliases
class GeneratorError(Exception):
    pass


class WorkerError(Exception):
    pass


class MissingJob(Exception):
    pass


# __________________________________________
class PoisonPill(object):
    "kill message"
    def __str__(self):return "PoisonPill"


# __________________________________________
def feed(q, g, m, verbose = False):
    """target of the jobgenerator, gets the jobs, measure the generation time,
       put the jobs into the input queue to the workers

       q = input queue
       g = job generator
       m = message queue"""
    nput, initime = 0, time.time()
    # --------
    jobid = 0
    while True:
        try:  
            start = time.time()
            job = g.next()
            gentime = (start, time.time())

        except StopIteration:
            if verbose:
                m.put(("InputQueue", time.time(), "put PoisonPill", None)) #***#

            # ending signal
            q.put(PoisonPill())
            return
        except Exception:
            # fatal error, the processing chain cannot survive to a generator error
            type, value, trace = sys.exc_info()
            message  = "JobGenerator could not generate job %d\n" % (jobid)
            message += "    " + "    ".join(traceback.format_exception(type, value, trace, limit=5))
            q.put(GeneratorError(message))
            return

        if verbose:
            m.put(("InputQueue", time.time(), "put job", jobid)) #***#
        q.put((jobid, job, gentime)); nput += 1 #count only jobs

        jobid += 1


# __________________________________________
class Job(object):
    """object to be yielded by the job generator
    arguments and keywordarguments will be passed to the worker in their respective workspaces
    """
    args, kwargs = (), {}
    def __init__(self, *args, **kwargs):
        self.args   = args
        self.kwargs = kwargs


# __________________________________________
class InputQueue(queues.Queue):
    """input queue joining the job generator and the workers"""
    def __init__(self, maxsize, generator, messagequeue):
        """
        :param maxsize: integer, maximum number of jobs that can fit in the inputqueue
        :param generator: generator, objects that can yield Job objects
        :param messagequeue: queues.Queue, queue where to leave verbose messages
        """
        queues.Queue.__init__(self, maxsize=maxsize)
        self.p = Process(target = feed, args = (self, generator, messagequeue, messagequeue is not None))
        self.p.start() 

    # ----------------------
    def close(self):
        self.p.join()        
        queues.Queue.close(self)


# __________________________________________
class BasicPrinter(Process):
    """Object to handle messages from the messagequeue
       the Basic Printer simply print messages on screen as they come from the messagequeue
    """
    def __init__(self, messagequeue):
        Process.__init__(self)
        self.messagequeue = messagequeue

    # ------------------------
    def communicate(self, message):
        print message

    # ------------------------
    def run(self):
        """custom section of any Process object"""
        while True:
            packet = self.messagequeue.get()
            if isinstance(packet, PoisonPill): break
            sender, tim, mess, jobid = packet
            message = "%s at %s : %s " % (sender + " " * (20 - len(sender)), str(ut(tim)), mess)
            if jobid is not None: message += str(jobid)
            print message #basic version : simply print on screen
        return

# __________________________________________
class NoPrinter(object):
    """the printer to use if the verbose mode is off"""
    def __init__(self, noqueue):
        self.pid = -1
        #messagequeue is None
    # ------------------------
    def start(self):
        return

    # ------------------------
    def join(self):
        return

    # ------------------------
    def terminate(self):
        return

    # ------------------------
    def communicate(self, message):
        print message


# __________________________________________
class ProcessPrinter(InteractiveStdOut):
    """the printer to use in process-oriented display
       all messages related to a given worker will appear on the same line on terminal
    """
    def __init__(self, messagequeue):
        InteractiveStdOut.__init__(self, maxlines = 1000, message_queue = messagequeue)

    # ------------------------
    def communicate(self, message):
        self.message_queue.put(("User", time.time(), message, None))

    # ------------------------
    def interpretor(self, tup):
        """how the message from the messagequeue should be handled"""
        if isinstance(tup, PoisonPill): 
            return -1, "exit iso" #send the exit signal to the printer

        sender, tim, mess, jobid = tup

        if   sender == "InputQueue": 
            line    = 0
            message = "%s at %s : %s " % (sender + " " * (20 - len(sender)), str(ut(tim)), mess)
            if jobid is not None: message += str(jobid)
        elif sender.split('-')[0] == "Worker": 
            line = int(sender.split('-')[-1])
            message = "%s at %s : %s " % (sender + " " * (20 - len(sender)), str(ut(tim)), mess)
            if jobid is not None: message += str(jobid)
        elif sender == "MessageQueue": 
            line = -1
            message = sender + mess
        elif sender == "User": 
            line = -1
            message = mess
        else: #raise Exception('message not understood')
            line = -1
            message = mess
        return line, message


# __________________________________________
class JobPrinter(InteractiveStdOut):
    """the printer to use in job-oriented display
       all messages related to a given job will appear on the same line on terminal
    """
    def __init__(self, messagequeue):
        InteractiveStdOut.__init__(self, maxlines = 100000, message_queue = messagequeue)

    # ------------------------
    def communicate(self, message):
        self.message_queue.put(("User", time.time(), message, None))

    # ------------------------
    def interpretor(self, tup):
        if isinstance(tup, PoisonPill): 
            return -1, "exit iso" #send the exit signal to the printer

        sender, tim, mess, jobid = tup
        if sender == "User": 
            line = -1
            message = mess
        elif jobid is None: 
            line = -1
            message = mess
        elif isinstance(jobid, int):
            line = jobid
            if sender == "InputQueue": message = "Job%d%s at %s  : %s" % (jobid, " " * (10 - len(str(jobid))), ut(tim), sender)
            elif "Worker"   in sender:
                if "got"    in mess: message = message = "Job%d%s at %s  : %s" % (jobid, " " * (10 - len(str(jobid))), ut(tim), sender)
                elif "put"  in mess: message = message = "Job%d%s at %s  : %s" % (jobid, " " * (10 - len(str(jobid))), ut(tim), "done")
                elif "fail" in mess: message = message = "Job%d%s at %s  : %s" % (jobid, " " * (10 - len(str(jobid))), ut(tim), "failed (see /tmp/multiproerrors.log)")
                else: line, message = -1, mess
            else: line, message = -1, mess
        else: line, message = -1, mess
        return line, message


# __________________________________________
class Target(object):
    """An object to host the target function or callable object provided by user
       this object will be stored into the workers (one copy for each of them)
       the jobs will be taken from the input queue (only way to enter the workspace of the workers)
       and results will be put in the outputqueue

    """
    def __init__(self, smth):
        """
        smth : function or callable object to be packed into the target object
               smth will receive arguments as jobs
               the function must be define as any other function, *args and **kwargs allowed
               e.g. >>def fun(x, y, *args, **kwargs): ...
               e.g. >>class Fun(object):
                    >>    def __call__(self, x, y, *args, **kwargs): ...
               the user may require the worker to be passed to the function by naming the first argument "worker"
               e.g. >>def fun(worker, x, y, *args, **kwargs): ...
               e.g. >>class Fun(object):
                    >>    def __call__(self, worker, x, y, *args, **kwargs): ...

        """
        if isinstance(smth, types.FunctionType):
            self.core        = smth
            self.parent      = None
            self.passparent  = False
            self.passworker  = inspect.getargspec(self.core).args[0] == "worker"
        elif isinstance(smth, types.ObjectType):
            self.core        = smth.__call__#target.__getattribute__("__call__")
            self.parent      = smth
            self.passparent  = True
            self.passworker  = inspect.getargspec(self.core).args[:2] == ["self", "worker"]
        else: raise NotImplementedError('')

    # --------------------
    def __call__(self, *args, **kwargs):
        if self.passparent: 
            #return self.core(self.parent, *args, **kwargs) #doesnt work (?)
            return self.parent(*args, **kwargs)
        else: return self.core(*args, **kwargs)


# __________________________________________
class Worker(Process):
    """The object that will run into a separate workspace, taking jobs from the iputqueue,
       passing them to the target object, and put results into the output queue.
       It may also put messages into the messagequeue
       the worker has some methods that may be useful for the target function and can enter it through the keyword argument "worker"
    """
    def __init__(self, target, inputqueue, outputqueue, messagequeue, raiseiferror=False, seed=None):
        """
        :param target: function or callable, the function to be attached to the worker
        :param inputqueue: queues.Queue, the queue where to find jobs
        :param outputqueue: queues.Queue, the queue where to put results
        :param messagequeue: queues.Queue or None, the queue where to put messages
        :param raiseiferror: bool, whether the whole process should be interrputed in case of failure
        :param seed: integer or None, the seed to attach to the worker for random applications
        """
        Process.__init__(self)
        self.inputqueue   = inputqueue
        self.outputqueue  = outputqueue
        self.messagequeue = messagequeue
        self.target       = target
        self.raiseiferror = raiseiferror
        self.seed         = seed
        self.verbose      = messagequeue is not None#not NoPrinter

        ##------ attach random functions to the worker
        if self.seed is None:
            #seedtmp    = self.pid + 10 * int((time.time() * 1.e4) % 10.)
            #randfuntmp = random.Random(seedtmp).random
            #self.seed = int(1000. * randfuntmp())
            #time.sleep(0.1)
            raise Exception('seed is None is no longer supported')
        self.rand_  = random.Random(self.seed).random
        self.tgauss = np.linspace(-10., 10., 100) 
        self.Fgauss = 0.5 * (1. + erf((self.tgauss - 0.) / (1. * np.sqrt(2.)))) #repart fun of the normal pdf
        self.Fgauss[0], self.Fgauss[-1] = 0., 1.0
        ##------

    #----
    def rand(self, N = 1):
        """return a random number based on a uniform pdf between 0 and 1
        if N is 1 : returns a float
        elif N > 1: returns an array"""
        if N == 1: return self.rand_()
        else: return np.array([self.rand_() for i in xrange(N)])

    #----
    def randn(self, N = 1):
        """same as rand, but following a gaussian pdf (mean 0, std 1)"""
        return np.interp(self.rand(N), xp = self.Fgauss, fp = self.tgauss)
        
    #----    
    def communicate(self, message):#!#
        """put message into the message queue"""
        self.messagequeue.put((self.name, time.time(), message, None)) #!# 

    #----
    def run(self):
        """gets jobs from the inputqueue and runs it until 
           it gets the ending signal
        """

        #----- for statistics
        ngot, nput, nfail  = 0, 0, 0
        inittime = time.time()
        #-----

        while True:
            packet = self.inputqueue.get()
            if isinstance(packet, PoisonPill):  #got the ending signal
                if self.verbose: self.messagequeue.put((self.name, time.time(), "got PoisonPill", None)) #***#

                self.inputqueue.put(packet)  #resend the ending signal for the other workers
                self.outputqueue.put(packet) #ending signal
                return
            elif isinstance(packet, GeneratorError):   #generator has failed
                self.inputqueue.put(PoisonPill())    #send the ending signal for the other workers
                self.outputqueue.put(packet) #transmit the GeneratorError
                return


            jobid, job, gentime = packet
            assert isinstance(job, Job)
            ngot += 1 #count only jobs, not signals or errors

            if self.verbose: self.messagequeue.put((self.name, time.time(), "got job", jobid)) #***#


            try: 
                start = time.time()
                if self.target.passworker: 
                    answer = self.target(self, *job.args, **job.kwargs) #pass self (i.e. the worker to self.target as first argument)
                else:               
                    answer = self.target(*job.args, **job.kwargs) #call the target function here!!!
                jobtime = (start, time.time())
            except:
                nfail += 1
                type, value, trace = sys.exc_info()
                message  = "Worker %s failed during job %d\n" % (self.name, jobid)
                message += "    " + "    ".join(traceback.format_exception(type, value, trace, limit=5))
                self.outputqueue.put(WorkerError(message))
                if self.raiseiferror: 
                    self.inputqueue.put(PoisonPill())    #send the ending signal for the other workers
                    return #stop the execution
                if self.verbose: self.messagequeue.put((self.name, time.time(), "failed", jobid)) #***#
                continue #ignore the error and continue getting tasks
            self.outputqueue.put((jobid, answer, gentime, jobtime))
            nput += 1 #count only jobs, not signals or errors

            if self.verbose: self.messagequeue.put((self.name, time.time(), "put job", jobid)) #***#
#----------------------
class WorkerStacker(Worker):
    """
    same as Worker, but do not put results into the output queue, unless the ending signal has been received
    use it for instance for stacking processes
    see tutorials
    """
    #----
    def run(self):
        """gets jobs from the inputqueue and runs it until 
           it gets the ending signal
        """

        #----- for statistics
        ngot, nput, nfail  = 0, 0, 0
        #!#inittime = time.time()
        jobids, answer = [],None #!#
        Tgen = 0.#!#
        Tpro = 0.#!#
        #-----

        while True:
            packet = self.inputqueue.get()
            if isinstance(packet, PoisonPill):  #got the ending signal
                if self.verbose: self.messagequeue.put((self.name, time.time(), "got PoisonPill", None)) #***#
                if answer is not None:#!#
                    self.outputqueue.put((jobids, answer, Tgen, Tpro))#!#

                self.inputqueue.put(packet)  #resend the ending signal for the other workers
                self.outputqueue.put(packet) #ending signal
                return
            elif isinstance(packet, GeneratorError):   #generator has failed
                if answer is not None:#!#
                    self.outputqueue.put((jobids, answer, Tgen, Tpro))#!#

                self.inputqueue.put(PoisonPill())    #send the ending signal for the other workers
                self.outputqueue.put(packet) #transmit the GeneratorError
                return


            jobid, job, gentime = packet
            assert isinstance(job, Job)
            ngot += 1 #count only jobs, not signals or errors
            if self.verbose: self.messagequeue.put((self.name, time.time(), "got job", jobid)) #***#


            try: 
                start = time.time()
                if self.target.passworker: 
                    answer = self.target(self, *job.args, **job.kwargs) #pass self (i.e. the worker) to self.target as first argument
                else:               
                    answer = self.target(*job.args, **job.kwargs) #call the target function here!!!
                jobtime = (start, time.time())
                Tgen += gentime[1] - gentime[0]#!#
                Tpro += jobtime[1] - jobtime[0]#!#
                jobids.append(jobid)#!#
            except:
                nfail += 1
                type, value, trace = sys.exc_info()
                message  = "Worker %s failed during job %d\n" % (self.name, jobid)
                message += "    " + "    ".join(traceback.format_exception(type, value, trace, limit=5))
                self.outputqueue.put(WorkerError(message))
                if self.raiseiferror: 
                    self.inputqueue.put(PoisonPill())    #send the ending signal for the other workers
                    return #stop the execution
                if self.verbose: self.messagequeue.put((self.name, time.time(), "failed", jobid)) #***#
                continue #ignore the error and continue getting tasks
            #!#self.outputqueue.put((jobid, answer, gentime, jobtime))
            #!#nput += 1 #count only jobs, not signals or errors
            if self.verbose: self.messagequeue.put((self.name, time.time(), "put job", jobid)) #***#
#__________________________________
class MapAsync(object):
    whichworker = Worker
    def __init__(self, funobj, generator, Nworkers=None, RaiseIfError=True, Taskset = None, Verbose=False, LowPriority=False):
        """Run a processing session, takes input from the job generator, and generates results when requested
        Arguments :
            funobj      = function/callble-object or Nworkers-list of functions/callble-objects,
                          the function(s) must be defined as usual, *args and **kwargs allowed
                          e.g. >>def fun(x, y, *args, **kwargs):
                               >>    t = kwargs['t']
                               >>    return x + y + t
                          e.g. >>class Fun(object):
                               >>    def __call__(self, x, y, *args, **kwargs):
                               >>        ...
                          the user may require the worker to be passed to the function
                          by naming the first argument "worker"
                          e.g. >>def fun(worker, x, y, *args, **kwargs):
                               >>    return x + y + kwargs['t'] + worker.randn()
                          e.g. >>class Fun(object): ...

            generator    = Job generator, must yield Job instances with arguments
                           and keyword arguments to be passed to the function
                           if used, do not include the worker argument in the jobs
                           e.g. (Job(x, y, t=1) for x, y in zip([1, 2], [2, 3]))
            Nworkers     = int, number of workers to use, None means use cpu_count
            RaiseIfError = bool, should we interrupt the whole process in case of failure
            Taskset      = string like "0-11" or None, the range of physical threads to be used (check with htop)
            LowPriority  = bool, run the processes with lower priority than default

        Use : always call it into a "with" statement, iterate over it to get results
            >>with MapAsync(myfunction, myjobgenerator, ...) as ma:
            >>    for jobid, results, generation_times, processing_times in ma:
            >>        #jobid = integer, index of the job as provided by the jobgenerator (starting with 0)
            >>        #        map async may return jobs in wrong order, see MapSync for ordered results
            >>        #results = whatever is returned by the function
            >>        #generation_times = tuple with two floats = starting and ending times at which the job was generated
            >>        #processing_times = tuple with two floats = starting and ending times at which the job was processed
            >>        do_something_with(results) #if too slow, then the parallelization may be inefficient

        See tutorials for detailed usage instructions

        """

        if Nworkers is None: Nworkers = cpu_count()
        # ----------- choose the message queue and printer
        self.Mes     = None
        self.printer = NoPrinter(self.Mes)
        if Verbose: #anything else than False
            self.Mes     = Queue(maxsize = 10000) #message queue  
            if Verbose is True:                 self.printer = BasicPrinter(self.Mes)
            elif Verbose.lower() == "job":      self.printer = JobPrinter(self.Mes)
            elif Verbose.lower() == "process":  self.printer = ProcessPrinter(self.Mes)
            elif Verbose.lower() == "basic":    self.printer = BasicPrinter(self.Mes)
            else:                               self.printer = BasicPrinter(self.Mes)

        # ----------- create the input and output queues
        self.In      = InputQueue(maxsize = Nworkers, generator = generator, messagequeue = self.Mes) #queue for inputs
        self.Out     = Queue(maxsize = Nworkers) #queue for outputs 
        # ne pas augmenter maxsize : il y a tjrs une meilleur solution a trouver!

        # ----------- determine if each worker will have a distinct target or not
        multifunobj = isinstance(funobj, list)
        if multifunobj:  assert len(funobj) == Nworkers

        #__________________________________
        self.rand    = random.Random().random
        self.raiseiferror = RaiseIfError
        self.taskset = Taskset
        self.lowpriority = LowPriority
        self.Nactive = Nworkers
        self.ppid    = os.getpid()
        self.verbose = self.Mes is not None
        #__________________________________
        self.workers = []
        seedstmp = np.random.rand(Nworkers) * 100000
        if not multifunobj: trgt = Target(funobj)
        for i in xrange(Nworkers):
            if multifunobj: trgt = Target(funobj[i])
            w = self.whichworker(\
                target       = trgt,
                inputqueue   = self.In,
                outputqueue  = self.Out, 
                messagequeue = self.Mes,
                raiseiferror = self.raiseiferror, 
                seed         = seedstmp[i]) #in case two mapasync a run at the same time
            w.name = "Worker-%04d" % (i + 1)
            self.workers.append(w)
    #______________________________________
    def __str__(self):
        s = '----------------------------\n'
        s += "Parent pid = %d\n" % self.ppid
        s += "    Generator pid = %d\n" % self.In.p.pid
        s += "    Printer   pid = %d\n" % self.printer.pid
        for w in self.workers:
            s += "    %s  pid = %d; seed = %d\n" % (w.name, w.pid, w.seed) 
        return s
    #______________________________________
    def __enter__(self):

        for w in self.workers: w.start()

        self.ppid = self.workers[0]._parent_pid
        self.pids = [w.pid for w in self.workers]
        self.pids.append(self.In.p.pid)
        self.settask()
        self.renice()

        self.pids.append(self.printer.pid)
        self.printer.start()


        return self
    #______________________________________
    def __exit__(self, type, value, trace):
        #case 1 : no error, all outputs have been extracted from Out => join the workers
        #case 2 : no error, not all outputs extracted from Out => terminate
        #case 3 : some errors => terminate
        if type is None and self.Nactive == 0:

            for w in self.workers: w.join() #this might be blocking
            self.printer.join()

            self.In.close()
            self.Out.close()
            if self.verbose:queues.Queue.close(self.Mes)

        else:
            #either an error has occured or the user leaves too soon
            #if self.verbose:print "killing workers and queues"
            
            queues.Queue.close(self.In)
            if self.verbose:queues.Queue.close(self.Mes)
            self.Out.close()
            self.printer.terminate()
            for w in self.workers: w.terminate()
            self.In.p.terminate()
        #------------

    #______________________________________
    def settask(self):
        if self.taskset is None: return
        if "-" in self.taskset:
            corestart, coreend = [int(x) for x in self.taskset.split('-')]
            assert coreend > corestart >= 0
            cmd = "taskset -pc %d-%d %%d" % (corestart, coreend)
            cmd = "\n".join([cmd % pid for pid in self.pids])
        else:
            corestart = coreend = int(self.taskset)
            assert coreend == corestart >= 0
            cmd = "taskset -pc %d %%d" % (corestart)
            cmd = "\n".join([cmd % pid for pid in self.pids])
        os.system(cmd)
    # ______________________________________
    def renice(self):
        if self.lowpriority:
            cmd = "renice -n 10 -p %d"
            cmd = "\n".join([cmd % pid for pid in self.pids])
            os.system(cmd)
        else:
            pass
    #______________________________________
    def __iter__(self): return self
    #______________________________________
    def communicate(self, *args, **kwargs):
        self.printer.communicate(*args, **kwargs)
    #______________________________________
    def next(self):
        if not self.Nactive: raise StopIteration
        while self.Nactive :
            packet = self.Out.get()
            if isinstance(packet, PoisonPill):
                self.Nactive -= 1
                continue
            elif isinstance(packet, GeneratorError):
                raise packet
            elif isinstance(packet, WorkerError) :
                with open('/tmp/multiproerrors.log', 'a') as fid:
                    fid.write(str(packet) + "\n")

                if self.raiseiferror: 
                    self.communicate(str(packet))
                    raise packet

                continue
            else:
                return packet#tuple jobid, answer, jobtime

        if self.verbose:
            self.Mes.put(("MessageQueue", time.time(), "got PoisonPill", None))
            self.Mes.put(PoisonPill())
        raise StopIteration
########################################### 
class waitqueue():
    """the object to reorder jobs once processed to be used by MapSync"""
    def __init__(self, generator, limit = 1e50): 
        """
        generator must return jobid, something
        the jobid list must be exaustive from 0 to N
        """
        self.l = []
        self.currentjob = 0
        self.generator  = generator
        self.limit      = limit
    #______________________________________
    def __len__(self): 
        return len(self.l)
    #______________________________________
    def append(self, jobid, packet):
        if len(self) >= self.limit: 
            raise Exception('the waitqueue was full')

        if not len(self.l): 
            self.l.append((jobid, packet))
            return
        i = 0
        while self.l[i][0] < jobid:
            i += 1
            if i == len(self.l): break
        self.l.insert(i, (jobid, packet))
    #______________________________________
    def pop(self, index):
        jobid, packet = self.l.pop(index)
        return packet #jobid, answer, ...
    #______________________________________
    def __iter__(self): return self
    #______________________________________
    def next(self):
        if len(self) and self.l[0][0] == self.currentjob: #first one is the right one
            self.currentjob += 1 
            return self.pop(0)

        while True:
            try:   
                packet = self.generator.next()#jobid, answer, ...
                jobid  = packet[0] #first item of the packet must be the jobid
            except StopIteration: break
        
            if jobid > self.currentjob: #came too soon, put it back
                self.append(jobid, packet)

            elif jobid == self.currentjob: #this is the one
                self.currentjob += 1
                return packet#jobid, answer, ...

        if len(self): 
            #may append if some processes have failed
            #print "warning : job %s never showed up" % self.currentjob
            self.currentjob += 1
            return self.currentjob - 1, MissingJob()
        raise StopIteration
#__________________________________________
class StackAsync(MapAsync):
    whichworker = WorkerStacker    
    """usage : 
    class UserStacker(object):
        "the target function designed for internal stack"
        def __init__(self, v = 0.):
            self.v = v #the current value of the stack
        def __call__(self, worker, x):
            "the way input data should be stacked with current value"
            start = time.time()
            while time.time() - start < 0.1: 0 + 0
            self.v += x
            return self, worker.name
        def __add__(self, other):
            "a method to merge stackers once back to the serial section"
            assert isinstance(other, UserStacker)
            return Stacker(self.v + other.v)

    def JobGen():
        for i in xrange(1000):
            yield Job(i)

    s = Stacker(0.) #the stacker to be duplicated in each worker
    with StackAsync(s, JobGen(), Verbose = False) as sa:
        S = None
        for jobids, (s, wname), Tgen, Tpro in sa: #receive partial stacks
            print "%s stacked %6d jobs in %.6fs, result %f" % (wname, len(jobids), Tgen + Tpro, s.v)
            if S is None: S = s
            else: S = S + s #merge all partial stacks using method __add__
    print "Final sum %f" % S.v
    """
#__________________________________________
class MapSync(MapAsync):
    def __init__(self, *args, **kwargs):
        """same as MapSync except that the job results are generated in the right order"""
        if 'RaiseIfError' in kwargs.keys() and not kwargs['RaiseIfError']:
            raise NotImplementedError("MapSync doesn't support RaiseIfError = False")
        MapAsync.__init__(self, *args, **kwargs)

    def __iter__(self):
        #jobs that come up too soon are kept in a waiting queue to preserve the input order
        return waitqueue(self)
#__________________________________________
class FakeWorker(object):
    """simulate a worker object"""
    def __init__(self):
        self.name  = "FakeWorker"
        self.rand  = np.random.rand
        self.randn = np.random.randn
        self.pid   = os.getpid()

#__________________________________________
class FakePrinter(object):
    """simulate a printer object"""
    def communicate(self, message):print message   
#__________________________________________
class FakeMapAsync(object):
    """
    use it istead of MapAsync to switch parallel computing off
    process is not parallelized and jobs are returned in the right order anyway
    """

    def __init__(self, funobj, generator, Nworkers=1, RaiseIfError=True, Taskset = None, Verbose=False, LowPriority=False):
        self.generator  = generator
        self.raiseiferror = RaiseIfError
        self.printer      = FakePrinter()
        self.jobid = -1
        self.target = Target(copy.deepcopy(funobj)) #copy funobj to avoid attribute modification as it is for MapAsync

    def __enter__(self): return self
    def __exit__(self, type, value, trace):return
    def __iter__(self): return self
    def communicate(self, message):
        self.printer.communicate(message)
    def next(self): 
        try : 
            start = time.time()
            job = self.generator.next()
            gentime = (start, time.time())
            self.jobid += 1
        except StopIteration: raise 
        except Exception as e: 
            raise e
    
        try:            
            start = time.time()
            if self.target.passworker: 
                answer = self.target(FakeWorker(), *job.args, **job.kwargs)
            else:                 
                answer = self.target(*job.args, **job.kwargs)
            jobtime = (start, time.time())
            return self.jobid, answer, gentime, jobtime
        except Exception as e:
            print e
            if self.raiseiferror: raise
#__________________________________________
class FakeMapSync(FakeMapAsync):
    """
    use it istead of MapSync to switch parallel computing off
    process is not parallelized and jobs are returned in the right order anyway
    """
    pass
#__________________________________________
def workfor(t):
    "fake working fot t seconds, processor activity will raise"
    start = time.time()
    while time.time() - start < t:
        tete_a_toto = 0. + 0.
#__________________________________________
if __name__ == "__main__":
    print "see tutorials for demo"


