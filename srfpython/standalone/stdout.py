import sys
import curses
import time
import numpy as np
import multiprocessing as mp

# --------------------------------------
def countdown(message, t):
    l = "%s %4.0fs" % (message, t)
    sys.stdout.write(l)
    sys.stdout.flush()

    start = time.time()
    while time.time() - start < t:
        sys.stdout.write("\b" * len(l))                

        l = "%s %4.0fs" % (message, t - time.time() + start)
        sys.stdout.write(l)
        sys.stdout.flush()
        time.sleep(0.5)


# --------------------------------------
class waitbar(object):
    def __str__(self):
        s = "%s %s %6.2f%% %s" % (self.title, self.bars(), 100. * self.purcent, self.remain)
        return s 
    def bars(self):
        nbars   = int(round(self.purcent * self.width))
        nspaces = self.width - nbars
        #return "%s" % (unichr(0x2588) * nbars) + "%s" % (" " * nspaces) #does not work of some systems
        return "%s" % ("|" * nbars) + "%s" % (" " * nspaces) #does not work of some systems

    # ____________________________________
    def __init__(self, title = "", width = 40, reevaluatespeed = 5.0):
        self.title   = title
        self.width   = width
        self.purcent = 0.
        self.lastpurcent = 0.
        self.start   = time.time()
        self.time     = self.start
        self.lasttime = self.start
        self.speed   = 0.
        self.remain  = "unkn"
        self.string  = self.__str__()
        self.reevaluatespeed = reevaluatespeed

        sys.stdout.write("%s" % self.string)
        sys.stdout.flush()

    # ____________________________________
    def refresh(self, purcent):
        if self.purcent < 0. or self.purcent > 1.: self.close()
        #if purcent and purcent - 1. and abs(self.purcent - purcent) < 0.005:return

        self.purcent, self.time = purcent, time.time()
        if self.reevaluatespeed and self.time - self.lasttime > self.reevaluatespeed:#reevaluate speed every X seconds
            self.speed = (self.purcent - self.lastpurcent) / (self.time - self.lasttime) 
            self.lastpurcent, self.lasttime = self.purcent, self.time
        elif self.start == self.lasttime:
            self.speed = (self.purcent - self.lastpurcent) / (self.time - self.lasttime) 

        if self.speed:
            tremain = ((1. - self.purcent) / self.speed)
            d = int(tremain / 24. / 3600.)
            h = int(tremain % (24. * 3600.) / 3600.)
            m = int(tremain % 3600. / 60.)
            s = int(tremain % 60.)
            self.remain = "%2ds" % s
            if m or h or d: self.remain = "%2dmn%s" % (m, self.remain)
            if h or d: self.remain = "%2dh%s" % (h, self.remain)
            if d: self.remain = "%3dd%s" % (d, self.remain)

        else:
            self.remain = "unkn"

        #self.remain = "(~%s remaining)" % self.remain
        self.remain = " %s" % self.remain
        self.remain = self.remain + " " * (25 - len(self.remain))
        sys.stdout.write("\b" * len(self.string))


        self.string  = self.__str__()
        sys.stdout.write("%s" % self.string)
        sys.stdout.flush()

    #____________________________________
    def close(self):
        self.purcent = 1.0-1.e-20
        self.refresh(self.purcent)
        sys.stdout.write("\n")


# --------------------------------------
class waitbarpipe(waitbar):
    #____________________________________
    def __str__(self):
        s = "%s %s %6.2f%% %s" % (self.title, self.bars(), 100. * self.purcent, self.remain)
        return s 
    def bars(self):
        nbars   = int(round(self.purcent * self.width))
        nspaces = self.width - nbars
        return "%s" % ("|" * nbars) + "%s" % ("." * nspaces)


# --------------------------------------
class multiline(object):
    def __init__(self, maxlines):
        self.maxlines = maxlines
    #____________________________________
    def __enter__(self):
        self.win = curses.initscr()
        curses.noecho()
        curses.cbreak()    
        self.win.keypad(True)
        self.win.scrollok(True)
        self.lines = ["" for i in xrange(self.maxlines)]
        self.line0 = 0 #first line printed
        self.lastcommunication = ""
        self.reset_termsize()
        return self

    #____________________________________
    def __exit__(self, type, value, trace):
        curses.echo()
        curses.nocbreak()
        self.win.keypad(0)
        self.win.scrollok(False)
        curses.endwin()

    #____________________________________
    def reset_termsize(self):
        self.nmax, self.mmax = self.win.getmaxyx()
        self.nmax -= 1
        self.mmax -= 1

    #____________________________________
    def refresh(self):
        self.win.refresh()   

    #____________________________________
    def communicate(self, message):
        self.lastcommunication = message
        self.reset_termsize()
        self.win.addstr(self.nmax, 0, message[:self.mmax])
        self.win.clrtoeol()
        self.refresh()

    #____________________________________
    def write(self, i, message, refresh = True):
        self.lines[i] = message

        while True:
            if   i <  self.line0: return
            elif i >  self.line0 + self.nmax: return
            try:
                self.win.addstr(i - self.line0, 0, message[:self.mmax])
                self.win.clrtoeol()
                break
            except curses.error:
                self.reset_termsize()
                continue           

        if refresh: self.refresh()   

    #____________________________________
    def move(self, line0):
        self.line0 = line0
        for i in xrange(self.line0, self.line0 + self.nmax + 1):
            self.write(i, self.lines[i], refresh = False)
        self.communicate(self.lastcommunication)
        self.refresh()
            
    #____________________________________
    def pause(self):
        self.communicate("pause")
        return self.win.getstr()


# --------------------------------------
class ExitedISO(Exception): pass
class InteractiveStdOut(mp.Process):
    def __init__(self, maxlines, message_queue = None):
        mp.Process.__init__(self)
        if message_queue is None: self.message_queue = mp.Queue()    #default queue
        else:                     self.message_queue = message_queue #user defined queue
        self.maxlines = maxlines

    def write(self, line, message):
        if not self.is_alive(): raise ExitedISO('')
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
        with  multiline(maxlines = self.maxlines) as ml:
            ml.win.nodelay(True)
            autofollow = False
            while True: 
                ml.reset_termsize()
                #------------------------
                try:
                    line, message = self.interpretor(self.message_queue.get_nowait())
                    if message.lower() == "exit iso": 
                        lines = ml.lines
                        break #ending signal from outside

                    if    line == -1:   ml.communicate(message)
                    else:               
                        ml.write(line, message)
                        if autofollow: 
                            if line - ml.line0 >= ml.nmax: 
                                ml.move(np.min([np.max([0, line - 2]), ml.maxlines - ml.nmax - 1]))

                except mp.queues.Empty: pass  
                except KeyboardInterrupt: raise
                except Exception as Detail:
                    ml.communicate("%s"% Detail)

                #------------------------
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
                elif ch in (curses.KEY_NPAGE, ord(' '), 0x06):   # Ctrl-F
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
                #elif ch == curses.KEY_LEFT: continue
                # cursor right
                #elif ch == curses.KEY_RIGHT: continue
                # toggle .dot-files
                elif ch == 0x08:# Ctrl-H
                    autofollow = not autofollow #toggle autofollow mode
                    ml.communicate('autofollow : %s' % str(autofollow))
                    
                # quit
                elif ch in (ord('q'), ord('Q')): 
                    lines = ml.lines
                    break#, curses.KEY_F10, 0x03): 
                else: continue
        #recall what was printed
        for l in lines: 
            if len(l): print l


# --------------------------------------
def MultiPrint(gen, maxlines = 10000):
    iso = InteractiveStdOut(maxlines = maxlines)

    iso.start()
    for line, message in gen:
        try:   iso.write(line, message)
        except ExitedISO: break

    iso.communicate('done, press q')
    while iso.is_alive(): time.sleep(1.)
    iso.join()


# ____________________________________
if __name__ == "__main__":
    if True:
        w = waitbar('test')
        for i in xrange(100):
            time.sleep(0.1)
            w.refresh(i / 100.)
        #w.refresh(1.)
        w.close()
    time.sleep(2.)
    if True:
        def gen():
            for i in xrange(100000):
                line = int(i / 1000) #chose a line number
                yield line, "ligne %d %f" % (line, np.random.rand())
        
        MultiPrint(gen(), 100) #prepare 120 lines for printing
