class bcolors:
    "source : http://stackoverflow.com/questions/287871/print-in-terminal-with-colors-using-python"
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def printpurple(*args):
    for l in args:
        print("%s%s%s" % (bcolors.HEADER, l, bcolors.ENDC), '\n')

def printblue(*args):
    for l in args:
        print("%s%s%s" % (bcolors.OKBLUE, l, bcolors.ENDC), '\n')

def printgreen(*args):
    for l in args:
        print("%s%s%s" % (bcolors.OKGREEN, l, bcolors.ENDC), '\n')

def printyellow(*args):
    for l in args:
        print("%s%s%s" % (bcolors.WARNING, l, bcolors.ENDC), '\n')

def printred(*args):
    for l in args:
        print("%s%s%s" % (bcolors.FAIL, l, bcolors.ENDC), '\n')

def printbold(*args):
    for l in args:
        print("%s%s%s" % (bcolors.BOLD, l, bcolors.ENDC), '\n')

def printunderline(*args):
    for l in args:
        print("%s%s%s" % (bcolors.UNDERLINE, l, bcolors.ENDC), '\n')
