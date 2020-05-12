from multiprocessing import queues, Process, Queue, Array, Value, sharedctypes, Lock, RLock
import sys, os
import traceback, inspect, random, time, types, copy, signal
import curses
import numpy as np
from math import erf as merf

global debug
debug = False
# if debug:
#     from obsmax4 import printblue, printgreen, printyellow, printred, printpurple

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

ML 24-07-2018
new feature : attach shared variables to the mapper so that it can be used by several processes, 
should minimize memory requirements

ML  27-07-2018
new feature : add the posibility to use generator targets 
avoids storing too much data in the worker workspace, see GenAsync

ML  30-07-2018
new feature : add a waitqueue for to reorder outputs from the GenAsync
attach it to GenSync mapper
move some components of obsmax4.tools.stdout here, reorganize the module
the module is now independent of obsmax4 except in debug mode

ML 14-09-2018
new feature : attach a lock object to the worker

"""


# ---------------- functions
def ut(t):
    s = time.ctime(t)
    return s.split()[3]


def workfor(t):
    "fake working fot t seconds, processor activity will raise (see htop)"
    start = time.time()
    while time.time() - start < t:
        tete_a_toto = 0. + 0.


def erf(x):
    return np.asarray([merf(w) for w in x], float)


# ---------------- errors and signals
class GeneratorError(Exception):
    """denote an error that occured inside the job generator"""
    pass


class WorkerError(Exception):
    """denote an error that occured inside a worker workspace"""
    pass


class MissingJob(Exception):
    """used in waitqueus to mentionned that some jobs never showed up"""
    pass


class PoisonPill(object):
    """used a a killin signal for communication between processes
    """

    def __str__(self):
        return "PoisonPill"


class ExitedISO(Exception):
    """used for advanced printing options"""
    pass


# ---------------- printers
class multiline(object):
    """prints stuff dynamically on several lines at the same time"""

    def __init__(self, maxlines):
        self.maxlines = maxlines

    # ____________________________________
    def __enter__(self):
        self.win = curses.initscr()
        curses.noecho()
        curses.cbreak()
        self.win.keypad(True)
        self.win.scrollok(True)
        self.lines = ["" for i in xrange(self.maxlines)]
        self.line0 = 0  # first line printed
        self.lastcommunication = ""
        self.reset_termsize()
        return self

    # ____________________________________
    def __exit__(self, type, value, trace):
        curses.echo()
        curses.nocbreak()
        self.win.keypad(0)
        self.win.scrollok(False)
        curses.endwin()

    # ____________________________________
    def reset_termsize(self):
        self.nmax, self.mmax = self.win.getmaxyx()
        self.nmax -= 1
        self.mmax -= 1

    # ____________________________________
    def refresh(self):
        self.win.refresh()

    # ____________________________________
    def communicate(self, message):
        self.lastcommunication = message
        self.reset_termsize()
        self.win.addstr(self.nmax, 0, message[:self.mmax])
        self.win.clrtoeol()
        self.refresh()

    # ____________________________________
    def write(self, i, message, refresh=True):
        self.lines[i] = message

        while True:
            if i < self.line0:
                return
            elif i > self.line0 + self.nmax:
                return
            try:
                self.win.addstr(i - self.line0, 0, message[:self.mmax])
                self.win.clrtoeol()
                break
            except curses.error:
                self.reset_termsize()
                continue

        if refresh: self.refresh()

        # ____________________________________

    def move(self, line0):
        self.line0 = line0
        for i in xrange(self.line0, self.line0 + self.nmax + 1):
            self.write(i, self.lines[i], refresh=False)
        self.communicate(self.lastcommunication)
        self.refresh()

    # ____________________________________
    def pause(self):
        self.communicate("pause")
        return self.win.getstr()


class InteractiveStdOut(Process):

    def __init__(self, maxlines, message_queue=None):
        Process.__init__(self)
        if message_queue is None:
            self.message_queue = Queue()  # default queue
        else:
            self.message_queue = message_queue  # user defined queue
        self.maxlines = maxlines

    def write(self, line, message):
        if not self.is_alive():
            raise ExitedISO('')
        self.message_queue.put((line, message))

    def communicate(self, message):
        """to be customized, put a message that must be understood by interpretor so that the output of interpretor will be : -1, message
        message "exit iso" forces the printer to leave (equivalent to pressing "q")
        """
        self.message_queue.put((-1, message))

    def interpretor(self, tup):
        """to be customized, tell me how to convert the message_queue outputs into a tuple like (line, message)"""
        line, message = tup
        return line, message

    def run(self):
        with multiline(maxlines=self.maxlines) as ml:
            ml.win.nodelay(True)
            autofollow = False
            while True:
                ml.reset_termsize()
                # ------------------------
                try:
                    line, message = self.interpretor(self.message_queue.get_nowait())
                    if message.lower() == "exit iso":
                        lines = ml.lines
                        break  # ending signal from outside

                    if line == -1:
                        ml.communicate(message)
                    else:
                        ml.write(line, message)
                        if autofollow:
                            if line - ml.line0 >= ml.nmax:
                                ml.move(np.min([np.max([0, line - 2]), ml.maxlines - ml.nmax - 1]))

                except queues.Empty:
                    pass
                except KeyboardInterrupt:
                    raise
                except Exception as Detail:
                    ml.communicate("%s" % Detail)

                # ------------------------
                ch = ml.win.getch()
                # cursor up
                if ch in (ord('k'), ord('K'), curses.KEY_UP):
                    ml.move(max([0, ml.line0 - 1]))
                    continue
                # cursor down
                elif ch in (ord('j'), ord('j'), curses.KEY_DOWN):
                    ml.move(min([ml.maxlines - ml.nmax - 1, ml.line0 + 1]))
                    continue
                    # page previous
                elif ch in (curses.KEY_PPAGE, curses.KEY_BACKSPACE, 0x02):
                    ml.move(max([0, ml.line0 - ml.nmax]))
                    continue
                # page next
                elif ch in (curses.KEY_NPAGE, ord(' '), 0x06):  # Ctrl-F
                    ml.move(min([ml.maxlines - ml.nmax - 1, ml.line0 + ml.nmax]))
                    continue
                # home
                elif ch in (curses.KEY_HOME, 0x01):
                    ml.move(0)
                    continue
                # end
                elif ch in (curses.KEY_END, 0x05):
                    ml.move(ml.maxlines - ml.nmax - 1)
                    continue
                # enter
                elif ch in (10, 13):
                    ml.move(min([ml.maxlines - ml.nmax - 1, ml.line0 + ml.nmax]))
                    continue
                # resize
                elif ch == curses.KEY_RESIZE:
                    ml.reset_termsize()
                    ml.move(ml.line0)
                    continue
                # cursor left
                # elif ch == curses.KEY_LEFT: continue
                # cursor right
                # elif ch == curses.KEY_RIGHT: continue
                # toggle .dot-files
                elif ch == 0x08:  # Ctrl-H
                    autofollow = not autofollow  # toggle autofollow mode
                    ml.communicate('autofollow : %s' % str(autofollow))

                # quit
                elif ch in (ord('q'), ord('Q')):
                    lines = ml.lines
                    break  # , curses.KEY_F10, 0x03):
                else:
                    continue
        # recall what was printed
        for l in lines:
            if len(l):
                print l


# def multiprint(gen, maxlines = 10000):
#     iso = InteractiveStdOut(maxlines = maxlines)
#
#     iso.start()
#     for line, message in gen:
#         try:
#             iso.write(line, message)
#         except ExitedISO:
#             break
#
#     iso.communicate('done, press q')
#     while iso.is_alive(): time.sleep(1.)
#     iso.join()


class NoPrinter(object):
    """use that printer to shut the processes up"""

    def __init__(self, noqueue):
        self.pid = -1
        # messagequeue is None

    def start(self):
        return

    def join(self):
        return

    def terminate(self):
        return

    def communicate(self, message):
        print message


class BasicPrinter(Process):
    """standard printing to stdout"""

    def __init__(self, messagequeue):
        Process.__init__(self)
        self.messagequeue = messagequeue

    def communicate(self, message):
        print message

    def run(self):
        while True:
            packet = self.messagequeue.get()
            if isinstance(packet, PoisonPill): break
            sender, tim, mess, jobid = packet
            message = "%s at %s : %s " % (sender + " " * (20 - len(sender)), str(ut(tim)), mess)
            if jobid is not None: message += str(jobid)
            print message
        return


class ProcessPrinter(InteractiveStdOut):
    def __init__(self, messagequeue):
        InteractiveStdOut.__init__(self, maxlines=1000, message_queue=messagequeue)

    def communicate(self, message):
        self.message_queue.put(("User", time.time(), message, None))

    def interpretor(self, tup):
        if isinstance(tup, PoisonPill):
            return -1, "exit iso"  # send the exit signal to the printer

        sender, tim, mess, jobid = tup

        if sender == "InputQueue":
            line = 0
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
        else:  # raise Exception('message not understood')
            line = -1
            message = mess
        return line, message


class JobPrinter(InteractiveStdOut):
    def __init__(self, messagequeue):
        InteractiveStdOut.__init__(self, maxlines=100000, message_queue=messagequeue)

    def communicate(self, message):
        self.message_queue.put(("User", time.time(), message, None))

    def interpretor(self, tup):
        if isinstance(tup, PoisonPill):
            return -1, "exit iso"  # send the exit signal to the printer

        sender, tim, mess, jobid = tup
        if sender == "User":
            line = -1
            message = mess
        elif jobid is None:
            line = -1
            message = mess
        elif isinstance(jobid, int):
            line = jobid
            if sender == "InputQueue":
                message = "Job%d%s at %s  : %s" % (jobid, " " * (10 - len(str(jobid))), ut(tim), sender)
            elif "Worker" in sender:
                if "got" in mess:
                    message = message = "Job%d%s at %s  : %s" % (jobid, " " * (10 - len(str(jobid))), ut(tim), sender)
                elif "put" in mess:
                    message = message = "Job%d%s at %s  : %s" % (jobid, " " * (10 - len(str(jobid))), ut(tim), "done")
                elif "fail" in mess:
                    message = message = "Job%d%s at %s  : %s" % (
                    jobid, " " * (10 - len(str(jobid))), ut(tim), "failed (see /tmp/multiproerrors.log)")
                else:
                    line, message = -1, mess
            else:
                line, message = -1, mess
        else:
            line, message = -1, mess
        return line, message


class FakePrinter(object):
    def communicate(self, message):
        print message


# ----------------
def feed(q, g, m, verbose=False):
    """ the target of a InputQueue, this function will feed the inputqueue with jobs
        from the job generator in the same time as the workers are getting jobs
    :param q: input queue
    :param g: job generator
    :param m: message queue"""
    nput, initime = 0, time.time()

    jobid = 0
    while True:
        try:
            start = time.time()
            job = g.next()
            gentime = (start, time.time())

        except StopIteration:
            if verbose: m.put(("InputQueue", time.time(), "put PoisonPill", None))  # ***#

            q.put(PoisonPill())  # ending signal
            return
        except:  # fatal error, the processing chain cannot survive to a generator error
            type, value, trace = sys.exc_info()
            message = "JobGenerator could not generate job %d\n" % (jobid)
            message += "    " + "    ".join(traceback.format_exception(type, value, trace, limit=5))
            q.put(GeneratorError(message))
            return

        if verbose: m.put(("InputQueue", time.time(), "put job", jobid))  # ***#
        q.put((jobid, job, gentime));
        nput += 1  # count only jobs

        jobid += 1


# ########################## Exchanging Queues
class InputQueue(queues.Queue):
    def __init__(self, maxsize, generator, messagequeue):
        queues.Queue.__init__(self, maxsize=maxsize)
        self.p = Process(target=feed, args=(self, generator, messagequeue, messagequeue is not None))
        self.p.start()

    def close(self):
        self.p.join()
        queues.Queue.close(self)


class WaitQueue_old:
    """receive processing outputs in a random order
    make them wait until the right ones shows up
    returns them in the correct order

    packet = received from the outputqueue
            tuple like (jobid, answer, gentime, jobtime)
            jobid is the number attributed to the job when it was generated,
            we use that number to re-order the packets
    """

    def __init__(self, generator, limit=1e50):
        raise Exception('obsolet')
        """
        generator must return jobid, something
        the jobid list must be exaustive from 0 to N
        """
        self.l = []  # a list with waiting packets
        self.currentjob = 0  # index of the currently expexted jobid
        self.generator = generator  # a generator of packets
        self.limit = limit  # max size

    def __len__(self):
        return len(self.l)

    def append(self, jobid, packet):
        if len(self) >= self.limit:
            raise Exception('the {} was full'.format(self.__class__.__name__))

        if not len(self.l):
            # first packet received
            self.l.append((jobid, packet))
            return

        # place (jobid, packet) at the right place in self.l
        i = 0
        while self.l[i][0] < jobid:
            i += 1
            if i == len(self.l):
                break
        self.l.insert(i, (jobid, packet))

    def pop(self, index):
        "extract the packet in self.l[index], remove it from self.l"
        jobid, packet = self.l.pop(index)
        return packet

    def __iter__(self):
        return self

    def next(self):
        if len(self) and self.l[0][0] == self.currentjob:
            # the first item in self.l is the expected packet,
            # remove it from self.l and return it
            # increment the expected packet number
            self.currentjob += 1
            return self.pop(0)

        while True:
            try:
                # get the next packet from the generator (i.e. the outputqueue)
                packet = self.generator.next()
                # the first item of the packet must be the jobid
                jobid = packet[0]
            except StopIteration:
                break

            if jobid > self.currentjob:
                # this packet came too soon, store it and do not return it yet
                self.append(jobid, packet)

            elif jobid == self.currentjob:
                # got the right packet, move to the net one
                self.currentjob += 1
                return packet

        if len(self):
            # may append if some processes have failed
            # print "warning : job %s never showed up" % self.currentjob
            self.currentjob += 1
            return self.currentjob - 1, MissingJob()
        raise StopIteration


class WaitQueue(object):
    """receive outputs from the output queue in a random order
        make them wait until the right ones shows up
        returns them in the correct order

        packet = received from the outputqueue
                tuple like (jobid, answer, gentime, jobtime)
                jobid is the number attributed to the job when it was generated,
                we use that number to re-order the packets
        """

    def __init__(self, generator, limit=1e50):
        """
        generator must return jobid, something
        the jobid list must be exaustive from 0 to N
        """
        self.jobids = []
        self.packets = []
        self.currentjob = 0  # index of the currently expexted jobid
        self.generator = generator  # a generator of packets
        self.limit = limit  # max size

    def __len__(self):
        return len(self.packets)

    def append(self, jobid, packet):
        if len(self) >= self.limit:
            raise Exception('the {} was full'.format(self.__class__.__name__))

        if not len(self.packets):
            # first packet received
            # self.l.append((jobid, packet))
            self.jobids.append(jobid)
            self.packets.append(packet)
            return

        # place (jobid, packet) at the right place in self.l
        i = np.searchsorted(self.jobids, jobid)
        self.jobids.insert(i, jobid)
        self.packets.insert(i, packet)

    def pop(self):
        "extract the packet in self.l[0], remove it from self.l"
        self.jobids.pop(0)
        packet = self.packets.pop(0)
        return packet

    def __iter__(self):
        return self

    def next(self):
        if len(self) and self.jobids[0] == self.currentjob:
            # the first item in self.l is the expected packet,
            # remove it from self.l and return it
            # increment the expected packet number
            self.currentjob += 1
            return self.pop()

        while True:
            try:
                # get the next packet from the generator (i.e. the outputqueue)
                packet = self.generator.next()
                # the first item of the packet must be the jobid
                jobid = packet[0]
            except StopIteration:
                break

            if jobid > self.currentjob:
                # this packet came too soon, store it and do not return it yet
                self.append(jobid, packet)

            elif jobid == self.currentjob:
                # got the right packet, move to the net one
                self.currentjob += 1
                return packet

        if len(self):
            # may append if some processes have failed
            # print "warning : job %s never showed up" % self.currentjob
            self.currentjob += 1
            return self.currentjob - 1, MissingJob()
        raise StopIteration


class WaitQueueGen(WaitQueue):
    """receive outputs from the output queue in a random order
        make them wait until the right ones shows up
        returns them in the correct order

        packet = received from the outputqueue
                tuple like ((jobid, nitem), answer, gentime, jobtime)
                jobid is the number attributed to the job when it was generated,
                nitem is the iteration number in job jobid, -1 means StopIteration
                we use that number to re-order the packets
        """

    def __init__(self, generator, limit=1e50):
        """
        generator must return jobid, something
        the jobid list must be exaustive from 0 to N

        example :
        self.jobids = [ 1,      4,      7]
        self.nitems = [[0, 3], [2, 5], [0, 1, -1]]
        self.packets = same shape as self.nitems


        """
        self.jobids = []
        self.nitems = []
        self.packets = []
        self.currentjob = 0  # index of the currently expexted jobid
        self.currentitem = 0
        self.generator = generator  # a generator of packets
        self.limit = limit  # max size

    def __len__(self):
        return len(self.jobids)

    def append(self, jobid, nitem, packet):

        if len(self) >= self.limit:
            raise Exception('the {} was full'.format(self.__class__.__name__))

        if not len(self.packets):
            # first packet received
            self.jobids.append(jobid)
            self.nitems.append([nitem])
            self.packets.append([packet])
            return

        # place (jobid, packet) at the right place in self.l
        i = np.searchsorted(self.jobids, jobid)
        if i < len(self.jobids) and self.jobids[i] == jobid:
            # jobid already in self.jobids
            # j = np.searchsorted(self.nitems[i], nitem)
            # self.nitems[i].insert(j, nitem)
            # self.packets[i].insert(j, packet)
            assert nitem == -1 or nitem == self.nitems[i][-1] + 1
            self.nitems[i].append(nitem)
            self.packets[i].append(packet)
        else:
            self.jobids.insert(i, jobid)
            self.nitems.insert(i, [nitem])
            self.packets.insert(i, [packet])

    def pop(self):
        "extract the packet in self.l[index], remove it from self.l"
        # print ">>>", self.jobids
        # print ">>>", self.nitems
        if len(self.nitems[0]) == 1 and self.nitems[0][0] == -1:
            # only one item with nitem==-1 remaining, remove this jobid from the store
            jobid = self.jobids.pop(0)  # e.g. 3
            nitem = self.nitems.pop(0)[0]  # e.g. [12]
            packet = self.packets.pop(0)[0]
        else:
            jobid = self.jobids[0]
            nitem = self.nitems[0].pop(0)
            packet = self.packets[0].pop(0)
        return packet

    def __iter__(self):
        return self

    def next(self):
        if len(self.jobids) and self.jobids[0] == self.currentjob and len(
                self.nitems[0]):  # and self.nitems[0] == self.currentitem:
            # the first item in self.l is the expected packet,
            # remove it from self.l and return it
            # increment the expected packet number
            if self.nitems[0][0] == -1:
                self.currentitem = 0
                self.currentjob += 1
            else:
                self.currentitem += 1

            packet = self.pop()
            if debug:
                printyellow('return from store', packet)
            return packet

        while True:
            try:
                # get the next packet from the generator (i.e. the outputqueue)
                packet = self.generator.next()
                if debug:
                    printgreen("$", packet)

                # the first item of the packet must be the jobid
                jobid, nitem = packet[0]
            except StopIteration:
                break

            if jobid > self.currentjob:  # and nitem > self.currentitem:
                # this packet came too soon, store it and do not return it yet
                if debug:
                    printpurple("append", jobid, nitem, packet)

                self.append(jobid, nitem, packet)

            elif jobid == self.currentjob:
                # got the right packet, move to the next one

                if nitem == -1:
                    if debug:
                        printred("???", self.jobids)
                        printred("???", self.nitems)
                    if len(self.nitems) and not len(self.nitems[0]):
                        # remove empty lists corresponding to jobid (done)
                        self.jobids.pop(0)
                        self.nitems.pop(0)
                        self.packets.pop(0)
                    self.currentjob += 1
                    self.currentitem = 0

                elif nitem == self.currentitem:
                    self.currentitem += 1

                if debug:
                    printblue("return directly", jobid, nitem, packet)
                return packet
                # print "append", jobid, nitem, packet
                # self.append(jobid, nitem, packet)

        if len(self):
            raise Exception('missing jobs')
            # may append if some processes have failed
            # print "warning : job %s never showed up" % self.currentjob
            self.currentjob += 1
            return self.currentjob - 1, MissingJob()
        raise StopIteration


# ###################################
class Job(object):
    """ an object to store job arguments and keywordarguments
    use it to pack jobs out of the job-generator

    def jobgenerator(N):
        for n in xrange(N):
            yield Job(n, twotimesn = 2 * n)
    """
    args, kwargs = (), {}

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs


class Target(object):
    def __init__(self, smth):
        """
        :param smth: a function or object with __call__ method
        """
        if isinstance(smth, types.FunctionType):
            self.core = smth
            self.parent = None
            self.passparent = False
            self.passworker = inspect.getargspec(self.core).args[0] == "worker"
        elif isinstance(smth, types.ObjectType):
            self.core = smth.__call__
            self.parent = smth
            self.passparent = True
            self.passworker = inspect.getargspec(self.core).args[:2] == ["self", "worker"]
        else:
            raise NotImplementedError('')

    def __call__(self, *args, **kwargs):
        if self.passparent:
            # return self.core(self.parent, *args, **kwargs) #doesnt work (?)
            return self.parent(*args, **kwargs)
        else:
            return self.core(*args, **kwargs)


# ################################### WORKERS
class Worker(Process):
    def __init__(self, target, inputqueue, outputqueue, messagequeue, raiseiferror=False, seed=None, parent=None, lock=None):
        """
        :param target: a Target object, the function or object to call inside the dedicated workspaces
        :param inputqueue: a InputQueue object, the queue that transmit jobs from the main workspace inside the dedicated workspaces
        :param outputqueue:
        :param messagequeue:
        :param raiseiferror:
        :param seed:
        :param parent: the parent object, may be MapAsync, MapSync, StackAsync, ...
        :param lock:
        """
        Process.__init__(self)
        self.inputqueue = inputqueue
        self.outputqueue = outputqueue
        self.messagequeue = messagequeue
        self.target = target
        self.raiseiferror = raiseiferror
        self.seed = seed
        self.verbose = messagequeue is not None  # not NoPrinter
        self.parent = parent
        self.is_locked = False
        self.lock = lock

        # ------ attach random functions to the worker
        if self.seed is None:
            # seedtmp    = self.pid + 10 * int((time.time() * 1.e4) % 10.)
            # randfuntmp = random.Random(seedtmp).random
            # self.seed = int(1000. * randfuntmp())
            # time.sleep(0.1)
            raise Exception('')
        self.rand_ = random.Random(self.seed).random
        self.tgauss = np.linspace(-10., 10., 100)
        self.Fgauss = 0.5 * (1. + erf((self.tgauss - 0.) / (1. * np.sqrt(2.))))  # repart fun of the normal pdf
        self.Fgauss[0], self.Fgauss[-1] = 0., 1.0
        # ------

    def acquire(self):
        if self.is_locked:
            raise Exception('{} is already locked'.format(self.name))
        self.is_locked = True
        self.lock.acquire()

    def release(self):
        if not self.is_locked:
            raise Exception('{} is not locked'.format(self.name))
        self.is_locked = False
        self.lock.release()

    def rand(self, N=1):
        if N == 1:
            return self.rand_()
        else:
            return np.array([self.rand_() for i in xrange(N)])

    def randn(self, N=1):
        return np.interp(self.rand(N), xp=self.Fgauss, fp=self.tgauss)

    def communicate(self, message):  # !#
        self.messagequeue.put((self.name, time.time(), message, None))  # !#

    def run(self):
        """gets jobs from the inputqueue and runs it until
           it gets the ending signal
        """

        # ----- for statistics
        ngot, nput, nfail = 0, 0, 0
        inittime = time.time()
        # -----

        while True:
            packet = self.inputqueue.get()
            if isinstance(packet, PoisonPill):  # got the ending signal
                if self.verbose: self.messagequeue.put((self.name, time.time(), "got PoisonPill", None))  # ***#

                self.inputqueue.put(packet)  # resend the ending signal for the other workers
                self.outputqueue.put(packet)  # ending signal
                return
            elif isinstance(packet, GeneratorError):  # generator has failed
                self.inputqueue.put(PoisonPill())  # send the ending signal for the other workers
                self.outputqueue.put(packet)  # transmit the GeneratorError
                return

            jobid, job, gentime = packet
            assert isinstance(job, Job)
            ngot += 1  # count only jobs, not signals or errors

            if self.verbose: self.messagequeue.put((self.name, time.time(), "got job", jobid))  # ***#

            try:
                start = time.time()
                if self.target.passworker:
                    answer = self.target(self, *job.args,
                                         **job.kwargs)  # pass self (i.e. the worker to self.target as first argument)
                else:
                    answer = self.target(*job.args, **job.kwargs)  # call the target function here!!!
                jobtime = (start, time.time())

            except:
                nfail += 1
                type, value, trace = sys.exc_info()
                message = "Worker %s failed during job %d\n" % (self.name, jobid)
                message += "    " + "    ".join(traceback.format_exception(type, value, trace, limit=5))

                self.outputqueue.put(WorkerError(message))

                if self.raiseiferror:
                    self.inputqueue.put(PoisonPill())  # send the ending signal for the other workers
                    return  # stop the execution

                if self.verbose:
                    self.messagequeue.put((self.name, time.time(), "failed", jobid))  # ***#
                continue  # ignore the error and continue getting tasks

            self.outputqueue.put((jobid, answer, gentime, jobtime))
            nput += 1  # count only jobs, not signals or errors

            if self.verbose:
                self.messagequeue.put((self.name, time.time(), "put job", jobid))  # ***#


class WorkerStacker(Worker):
    """
    same as Worker, but do not put results into the output queue, unless the ending signal has been received
    use it for istance for stacking processes
    see tutorials
    """

    def run(self):
        """gets jobs from the inputqueue and runs it until
           it gets the ending signal
        """

        # ----- for statistics
        ngot, nput, nfail = 0, 0, 0
        # !#inittime = time.time()
        jobids, answer = [], None  # !#
        Tgen = 0.  # !#
        Tpro = 0.  # !#
        # -----

        while True:
            packet = self.inputqueue.get()
            if isinstance(packet, PoisonPill):  # got the ending signal
                if self.verbose: self.messagequeue.put((self.name, time.time(), "got PoisonPill", None))  # ***#
                if answer is not None:  # !#
                    self.outputqueue.put((jobids, answer, Tgen, Tpro))  # !#

                self.inputqueue.put(packet)  # resend the ending signal for the other workers
                self.outputqueue.put(packet)  # ending signal
                return
            elif isinstance(packet, GeneratorError):  # generator has failed
                if answer is not None:  # !#
                    self.outputqueue.put((jobids, answer, Tgen, Tpro))  # !#

                self.inputqueue.put(PoisonPill())  # send the ending signal for the other workers
                self.outputqueue.put(packet)  # transmit the GeneratorError
                return

            jobid, job, gentime = packet
            assert isinstance(job, Job)
            ngot += 1  # count only jobs, not signals or errors
            if self.verbose: self.messagequeue.put((self.name, time.time(), "got job", jobid))  # ***#

            try:
                start = time.time()
                if self.target.passworker:
                    answer = self.target(self, *job.args,
                                         **job.kwargs)  # pass self (i.e. the worker) to self.target as first argument
                else:
                    answer = self.target(*job.args, **job.kwargs)  # call the target function here!!!
                jobtime = (start, time.time())
                Tgen += gentime[1] - gentime[0]  # !#
                Tpro += jobtime[1] - jobtime[0]  # !#
                jobids.append(jobid)  # !#
            except:
                nfail += 1
                type, value, trace = sys.exc_info()
                message = "Worker %s failed during job %d\n" % (self.name, jobid)
                message += "    " + "    ".join(traceback.format_exception(type, value, trace, limit=5))
                self.outputqueue.put(WorkerError(message))
                if self.raiseiferror:
                    self.inputqueue.put(PoisonPill())  # send the ending signal for the other workers
                    return  # stop the execution
                if self.verbose: self.messagequeue.put((self.name, time.time(), "failed", jobid))  # ***#
                continue  # ignore the error and continue getting tasks
            # !#self.outputqueue.put((jobid, answer, gentime, jobtime))
            # !#nput += 1 #count only jobs, not signals or errors
            if self.verbose: self.messagequeue.put((self.name, time.time(), "put job", jobid))  # ***#


class WorkerGenerator(Worker):
    """
    same as Worker, but do not put results into the output queue, unless the ending signal has been received
    use it for istance for stacking processes
    see tutorials
    """

    def run(self):
        """gets jobs from the inputqueue and runs it until
           it gets the ending signal
        """

        # ----- for statistics
        ngot, nput, nfail = 0, 0, 0
        inittime = time.time()
        # -----

        while True:
            packet = self.inputqueue.get()
            if isinstance(packet, PoisonPill):  # got the ending signal
                if self.verbose:
                    self.messagequeue.put((self.name, time.time(), "got PoisonPill", None))  # ***#

                self.inputqueue.put(packet)  # resend the ending signal for the other workers
                self.outputqueue.put(packet)  # ending signal
                return
            elif isinstance(packet, GeneratorError):  # generator has failed
                self.inputqueue.put(PoisonPill())  # send the ending signal for the other workers
                self.outputqueue.put(packet)  # transmit the GeneratorError
                return

            jobid, job, gentime = packet
            assert isinstance(job, Job)
            ngot += 1  # count only jobs, not signals or errors

            if self.verbose:
                self.messagequeue.put((self.name, time.time(), "got job", jobid))  # ***#

            try:
                start = time.time()
                if self.target.passworker:
                    answer_generator = self.target(self, *job.args,
                                                   **job.kwargs)  # pass self (i.e. the worker to self.target as first argument)
                else:
                    answer_generator = self.target(*job.args, **job.kwargs)  # call the target function here!!!

                for nanswer, answer in enumerate(answer_generator):
                    self.outputqueue.put(((jobid, nanswer), answer, (-1., -1.), (-1., -1.)))
                jobtime = (start, time.time())
            except Exception:
                nfail += 1
                type, value, trace = sys.exc_info()
                message = "Worker %s failed during job %d\n" % (self.name, jobid)
                message += "    " + "    ".join(traceback.format_exception(type, value, trace, limit=5))
                self.outputqueue.put(WorkerError(message))
                if self.raiseiferror:
                    self.inputqueue.put(PoisonPill())  # send the ending signal for the other workers
                    return  # stop the execution
                if self.verbose:
                    self.messagequeue.put((self.name, time.time(), "failed", jobid))  # ***#
                continue  # ignore the error and continue getting tasks

            # put one more item in the queue to notify that the iteration is over
            # use nitem = -1
            # use answer = StopIteration instance
            self.outputqueue.put(((jobid, -1), StopIteration('job {} done'.format(jobid)), gentime, jobtime))
            nput += 1  # count only jobs, not signals or errors

            if self.verbose:
                self.messagequeue.put((self.name, time.time(), "put job", jobid))  # ***#


class FakeWorker(object):
    def __init__(self):
        self.name = "FakeWorker"
        self.rand = np.random.rand
        self.randn = np.random.randn
        self.pid = os.getpid()

    def acquire(self):
        pass

    def release(self):
        pass


class FakeLock(object):
    def acquire(self):
        pass

    def release(self):
        pass


# ################################### MAPPERS
class MapAsync(object):
    whichworker = Worker

    def __init__(self, funobj, generator, Nworkers=12, RaiseIfError=True, Taskset=None,
                 Verbose=False, LowPriority=False,
                 SharedVariables=None):
        """
        funobj    = function/callble-object or Nworkers-list of functions/callble-objects, may include the "worker" argument as first argument
        generator = Job generator, must yield Job instances with arguments and keyword arguments to be passed
        Nworkers  = int, number of workers to use
        ...
        """
        if Nworkers is None:
            from multiprocessing import cpu_count
            Nworkers = cpu_count()
        # ----------- choose the message queue and printer
        self.Mes = None
        self.printer = NoPrinter(self.Mes)
        if Verbose:  # anything else than False
            self.Mes = Queue(maxsize=10000)  # message queue
            if Verbose is True:
                self.printer = BasicPrinter(self.Mes)
            elif Verbose.lower() == "job":
                self.printer = JobPrinter(self.Mes)
            elif Verbose.lower() == "process":
                self.printer = ProcessPrinter(self.Mes)
            elif Verbose.lower() == "basic":
                self.printer = BasicPrinter(self.Mes)
            else:
                self.printer = BasicPrinter(self.Mes)

        # ----------- create the input and output queues
        self.In = InputQueue(maxsize=Nworkers, generator=generator, messagequeue=self.Mes)  # queue for inputs
        self.Out = Queue(maxsize=Nworkers)  # queue for outputs
        # ne pas augmenter maxsize : il y a tjrs une meilleur solution a trouver!

        # ----------- determine if each worker will have a distinct target or not
        multifunobj = isinstance(funobj, list)
        if multifunobj:
            assert len(funobj) == Nworkers

        self.rand = random.Random().random
        self.raiseiferror = RaiseIfError
        self.taskset = Taskset
        self.lowpriority = LowPriority
        self.Nactive = Nworkers
        self.ppid = os.getpid()
        self.verbose = self.Mes is not None
        self.lock = RLock()  # one locker shared between all workers

        # attach shared variables to self
        if SharedVariables is not None:
            for key, val in SharedVariables.items():
                assert isinstance(val, sharedctypes.SynchronizedArray) or \
                       isinstance(val, sharedctypes.Synchronized)
                self.__setattr__(key, val)

        self.workers = []
        seedstmp = np.random.rand(Nworkers) * 100000
        if not multifunobj: trgt = Target(funobj)
        for i in range(Nworkers):
            if multifunobj: trgt = Target(funobj[i])
            w = self.whichworker( \
                target=trgt,
                inputqueue=self.In,
                outputqueue=self.Out,
                messagequeue=self.Mes,
                raiseiferror=self.raiseiferror,
                seed=seedstmp[i],  # in case two mapasync run at the same time
                parent=self,
                lock=self.lock)
            w.name = "Worker-%04d" % (i + 1)
            self.workers.append(w)

    def __str__(self):
        s = '----------------------------\n'
        s += "Parent pid = %d\n" % self.ppid
        s += "    Generator pid = %d\n" % self.In.p.pid
        s += "    Printer   pid = %d\n" % self.printer.pid
        for w in self.workers:
            s += "    %s  pid = %d; seed = %d\n" % (w.name, w.pid, w.seed)
        return s

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

    def __exit__(self, type, value, trace):
        # case 1 : no error, all outputs have been extracted from Out => join the workers
        # case 2 : no error, not all outputs extracted from Out => terminate
        # case 3 : some errors => terminate
        if type is None and self.Nactive == 0:

            for w in self.workers: w.join()  # this might be blocking
            self.printer.join()

            self.In.close()
            self.Out.close()
            if self.verbose: queues.Queue.close(self.Mes)

        else:
            # either an error has occured or the user leaves too soon
            # if self.verbose:print "killing workers and queues"

            queues.Queue.close(self.In)
            if self.verbose: queues.Queue.close(self.Mes)
            self.Out.close()
            self.printer.terminate()
            for w in self.workers: w.terminate()
            self.In.p.terminate()

    def settask(self):
        if self.taskset is None: return
        if "-" in self.taskset:
            corestart, coreend = [int(x) for x in self.taskset.split('-')]
            assert coreend > corestart >= 0
            cmd = "taskset -pca %d-%d %%d" % (corestart, coreend)
            cmd = "\n".join([cmd % pid for pid in self.pids])
        else:
            corestart = coreend = int(self.taskset)
            assert coreend == corestart >= 0
            cmd = "taskset -pca %d %%d" % (corestart)
            cmd = "\n".join([cmd % pid for pid in self.pids])
        os.system(cmd)

    def renice(self):
        if self.lowpriority:
            # cmd = "renice -n 10 -p %d"
            # cmd = "\n".join([cmd % pid for pid in self.pids])
            cmd = "renice -n 10 -g %d" % self.ppid
            os.system(cmd)
        else:
            pass

    def __iter__(self):
        return self

    def communicate(self, *args, **kwargs):
        self.printer.communicate(*args, **kwargs)

    def next(self):
        if not self.Nactive:
            raise StopIteration
        while self.Nactive:
            packet = self.Out.get()
            if isinstance(packet, PoisonPill):
                self.Nactive -= 1
                continue
            elif isinstance(packet, GeneratorError):
                raise packet
            elif isinstance(packet, WorkerError):
                with open('multiproerrors.log', 'a') as fid:
                    fid.write(str(packet) + "\n")

                self.communicate(str(packet))
                if self.raiseiferror:
                    raise packet

                continue
            else:
                return packet  # tuple jobid, answer, jobtime

        if self.verbose:
            self.Mes.put(("MessageQueue", time.time(), "got PoisonPill", None))
            self.Mes.put(PoisonPill())
        raise StopIteration


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


class GenAsync(MapAsync):
    whichworker = WorkerGenerator
    """usage : 

    def JobGen():
        for jobid in xrange(10):
            yield Job(N=10)

    def Target(N):
        # must be a generator, no check
        for nitem in xrange(N):
            answer = np.cos(nitem)
            yield answer

    with GenAsync(Target, JobGen()) as ga:
        for (jobid, nitem), answer, _, _ in ga:
            # jobid is the number of the job
            # nitem is the number of the item yielded by the target function that processed job jobid
            print jobid, nitem, answer

    """


class MapSync(MapAsync):
    whichwaitqueue = WaitQueue

    def __init__(self, *args, **kwargs):
        if 'RaiseIfError' in kwargs.keys() and not kwargs['RaiseIfError']:
            raise NotImplementedError("MapSync doesn't support RaiseIfError = False")
        MapAsync.__init__(self, *args, **kwargs)

    def __iter__(self):
        # jobs that come up too soon are kept in a waiting queue to preserve the input order
        return self.whichwaitqueue(self)


class GenSync(MapSync):
    whichworker = WorkerGenerator
    whichwaitqueue = WaitQueueGen  # : need a waitqueue that can handle (jobid, nitem) instead of jobid


class FakeMapAsync(object):
    """
    use it istead of MapAsync to switch parallel computing off
    process is not parallelized and jobs are returned in the right order anyway
    """

    def __init__(self, funobj, generator, Nworkers=1, RaiseIfError=True, Taskset=None, Verbose=False,
                 LowPriority=False):
        self.generator = generator
        self.raiseiferror = RaiseIfError
        self.printer = FakePrinter()
        self.jobid = -1
        self.target = Target(copy.deepcopy(funobj))  # copy funobj to avoid attribute modification as it is for MapAsync

    def __enter__(self):
        return self

    def __exit__(self, type, value, trace):
        return

    def __iter__(self):
        return self

    def communicate(self, message):
        self.printer.communicate(message)

    def next(self):
        try:
            start = time.time()
            job = self.generator.next()
            gentime = (start, time.time())
            self.jobid += 1
        except StopIteration:
            raise
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


class FakeMapSync(FakeMapAsync):
    """
    use it istead of MapSync to switch parallel computing off
    process is not parallelized and jobs are returned in the right order anyway
    """
    pass


if __name__ == "__main__":
    from obsmax4.graphictools.gutils import *
    from obsmax4.graphictools.cmaps import *
    from obsmax4.statstools.pdf import PDF, ccl
    import numpy as np
    import os

    Njob = 96
    Nworkers = 12
    # ----------------------
    if False:  # 1) use a custom function
        def fun(worker, n, i):
            """Define your own target function, use arguments or keyword arguments
            that will be passed from the job generator and return outputs that will be
            pushed to the output queue.
            Using argument "worker" as first argument make the worker be passed to target function
            for thread safe applications (i.e. random numbers)"""
            workfor(1.0 + 2.0 * worker.rand())
            # if i == 3: raise Exception('')
            return worker.pid, worker.rand(n)
    else:  # 2) use a custom object with dedicated attributes to be stored into the target
        class Obj(object):
            def __init__(self, cst=0):
                """put here the target attributes
                   remember that the process runs into a separate workspace,
                   modifications done by the target will be lost after closing the mapper!!
                   this should only be used to avoid too large data exchanges by
                   packing into the target the variables that will be used by all jobs"""
                self.cst = cst

            def __call__(self, worker, n, i):
                """modifying the "self" attributes will only have an effect inside the dedicated workspace!!!
                """
                np.median(np.random.randn(800000))
                workfor(1.01 + 2.02 * worker.rand())  # workfor(1.0 + 2.0 * worker.rand())
                # if i == 3: raise Exception('')
                return worker.pid, worker.rand(n) + self.cst * 5.


        # fun = Obj()
        fun = [Obj(cst=i) for i in xrange(Nworkers)]  # one specific callable per worker


    # ----------------------
    def gen(n):
        """
        create the list of arguments and keyword arguments
        use the Job class to pack the arguments in a convenient package that will be passed to the
        target as input arguments
        use the yield command instead of return, so the jobs will be generated
        when required by the workers
        """
        for i in xrange(n):
            workfor(0.1)
            # if i == 3: raise Exception('')
            yield Job(n=10, i=i)


    # ----------------------
    lgd, pids = [], []
    fig1 = plt.figure()
    fig2 = plt.figure()
    fig3 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(211)
    ax3 = fig2.add_subplot(212, sharex=ax2)

    # ----------------------
    with MapSync(fun, gen(Njob), Nworkers=Nworkers, Taskset='0-11',
                 Verbose=True) as ma:  # Verbose = "Process" Verbose = "Job"
        """call the mapper with the command with (unindent means closing the processors)
        use the MapAsync      mapper to get the results as they come (un ordered)
        --- --- MapSync       ------ -- --- --- ------- in correct order
        --- --- FakeMapAsync  ------ -- switch parallel computing off
        """
        print ma
        randomseries = []
        for jobid, (workerid, x), (genstart, genend), (jobstart, jobend) in ma:
            """iterate over the mapper to get the results"""
            # print "jobid%d adv%d" % (jobid, int(round((jobend - jobstart) / (genend - genstart))) + 1) #advised worker number
            ma.printer.communicate("jobid-%d adviced-worker-number-%d" % (
            jobid, int(round((jobend - jobstart) / (genend - genstart))) + 1))

            clr = np.random.rand(3)  # value2color(jobid / float(Njob), cmap = linecmap(Njob))
            ax1.plot(x, color=clr)
            ax2.plot([genstart, genend], [jobid, jobid], "k", linewidth=3)
            ax2.plot([jobstart, jobend], [jobid, jobid], color=clr, linewidth=3)
            ax2.text(jobend, jobid, "pid %d" % workerid)
            ax3.plot([jobstart, jobend], [workerid, workerid], color=clr, linewidth=3)

            pids.append(workerid)

            lgd.append("pid %d" % workerid)
            randomseries.append(x)

    # ----------------------
    ax1.legend(lgd)
    timetick(ax2, "x")
    ax2.grid()

    timetick(ax3, "x")
    pids = np.sort(np.unique(pids))
    ax3.set_yticks(pids)
    ax3.set_yticklabels(pids)
    ax3.grid()

    ax1.figure.show()
    ax2.figure.show()

    #    C = np.zeros((len(randomseries), len(randomseries)))
    #    for i in xrange(len(randomseries)):
    #        for j in xrange(i, len(randomseries)):
    #            C[i, j] = C[j, i] = ccl(randomseries[i], randomseries[j])

    #    plt.figure()
    #    gca().pcolormesh1(C, vmin = -1.0, vmax = 1.0)
    #    gcf().show()

    pause()

