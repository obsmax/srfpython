from __future__ import print_function
import time
import signal
import imp
import os
import sys
import numpy as np
from numpy import log10, inf, nan #needed by eval


class TimeOutError(Exception):
    pass


class Timeout():

    def __init__(self, sec):
        self.sec = sec

    def __enter__(self):
        signal.signal(signal.SIGALRM, self.raise_timeout)
        signal.alarm(self.sec)

    def __exit__(self, *args):
        signal.alarm(0)

    def raise_timeout(self, *args):
        raise TimeOutError()


class Timer(object):

    def __init__(self, title):
        self.title = title

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args, **kwargs):
        print ("elapsed time %s : %fs" % (self.title, time.time() - self.start))


def firstfalse(I):
    if not I[0] : return 0
    II = I[1:] != I[:-1]
    r = np.arange(len(I) - 1)
    return r[II][0] + 1


def munique(*Xs):
    assert np.all(len(X) == len(Xs[0]) for X in Xs)
    L = []
    for tup in zip(*Xs):
        if tup not in L:
            L.append(tup)
    return tuple([np.array(w) for w in zip(*L)])


def minmax(X):
    if hasattr(X, "min"): #numpy arrays
        return X.min(), X.max()
    else:
        return min(X), max(X)


def discrete_time_primitive(x, u, area=False):
    if len(x) - len(u):
        raise Exception('shape mismatch')
    if (x[1:] < x[:-1]).any():
        raise Exception('x must be sorted asc')

    v = np.zeros(len(u), u.dtype)
    v[1:] = 0.5 * (u[1:] + u[:-1]) * (x[1:] - x[:-1])
    if area:
        return np.sum(v)  # area below the curve
    else:
        return v.cumsum()  # primitive function to be associated to x


def cosTaperwidth(data, sampling_rate, width):
    tap = np.ones_like(data)
    Nwidth = int(np.round(width * sampling_rate))
    if not Nwidth : return tap

    t = np.arange(len(data)) / sampling_rate
    ttap = 0.5 * (np.sin(np.pi * t[:Nwidth] / float(width) / 1.0 + np.pi / 2.0) + 1.)
    tap[:Nwidth] *= ttap[::-1]
    tap[-Nwidth:] *= ttap
    return tap


# -------------------------------------
def string2func(s):
    "converts string into callable function using temporary python file"
    pyfilename = "/tmp/%s.py" % randstr(10)
    with open(pyfilename, 'w') as fid: fid.write(s)
    funcname = s.split('def ')[-1].split('(')[0].strip()
    func = getattr(imp.load_source(pyfilename.split('/')[-1].split('.py')[0], pyfilename), funcname)
    os.remove(pyfilename)
    os.remove(pyfilename + "c")
    return func


# -------------------------------------
def tostr(l, fmt):
    return " ".join(fmt % v for v in l)


# -------------------------------------
def randstr(n):
    """generate random strings with letters and numbers, not thread safe"""
    chrs   = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
    indexs = np.floor(np.random.rand(n) * len(chrs))
    return "".join([chrs[int(i)] for i in indexs])


# -------------------------------------
def string2func(s):
    "converts string into callable function using temporary python file"
    pyfilename = "/tmp/%s.py" % randstr(10)
    with open(pyfilename, 'w') as fid: fid.write(s)
    funcname = s.split('def ')[-1].split('(')[0].strip()
    func = getattr(imp.load_source(pyfilename.split('/')[-1].split('.py')[0], pyfilename), funcname)
    os.remove(pyfilename)
    os.remove(pyfilename + "c")
    return func


# --------------------------------------
def readargv():
    "read sys.argv and store results into a dictionary"
    def isnumeric(a):
        try:
            float(a)
            return True
        except ValueError:
            return False

    l = sys.argv[1:] #not array!
    for n, arg in enumerate(l):
        if n and arg[0] == "-" and not isnumeric(arg):
            l[n] = "__qwerty__%s" % arg

    l = "__azerty__".join(l)
    l = l.split('__azerty____qwerty__')
    l = [w.split('__azerty__') for w in l]

    D = {}#structure({})
    keyorder = []
    for ll in l:
        if ll[0][0] == "-":
            key = ll[0].strip('-')
            if key in D.keys(): raise Exception('argument repeated %s' % key)
            D[key] = ll[1:]
            keyorder.append(key)
        elif not hasattr(D, 'remain'):
            D['remain'] = ll
            keyorder.append('remain')
        else:
            raise Exception('this should not append')
    for k, v in D.items():
        for n, vv in enumerate(v):
            if isnumeric(vv): D[k][n] = eval(vv)
        if len(v) > 1:        D[k] = np.array(v)
        elif len(v) == 1:
            D[k] = [v[0]] #not [v[0]]
        else:
            D[k] = None
    D['_keyorder'] = keyorder
    return D

# --------------------------------------
def readargv1():
    """"read sys.argv and store results into a dictionary
    assume the argument list looks like

    1 2 3 \
    -option0 1 \
    --plugin1 a 1.0 -2.0  \
        -option1 ./*py \
        -option2 1 2. kjhkjh \
    --plugin2


    => returns :
        {'main': [1, 2, 3],
         'option0': [1],
         'plugin1': {'main': ['a', 1.0, -2.0],
                    'option1': ['./__init__.py', './utils.py'],
                    'option2': [1, 2.0, 'kjhkjh'],
                    '_keyorder': ['option1', 'option2']
                    }
         'plugin2': {'main': [],
                    '_keyorder': []
                    },
         '_keyorder': ['option0', 'plugin1', 'plugin2']
         }


    """

    # -----
    def isnumeric(a):
        try:
            float(a)
            return True
        except ValueError:
            return False

    # -----
    def evalifpossible(a):
        if isnumeric(a):
            try:
                return eval(a)
            except Exception:
                return a
        return a

    # -----
    def isplugin(arg):
        return arg.startswith('--')

    def isoption(arg):
        return not isplugin(arg) \
               and arg.startswith('-') \
               and not isnumeric(arg)

    def readplugin(l, parent=None):
        """recursive function"""
        plugin = {"main": [], "_keyorder": []}
        if not len(l):
            return plugin, l

        while len(l):
            arg = l[0]

            if isplugin(arg):
                if parent is None: #highest level
                    pluginname = arg#.split('-')[-1]
                    plugin[pluginname], l = readplugin(l[1:], parent=plugin)
                    plugin['_keyorder'].append(pluginname)
                else: #already inside a plugin, just leave the current one
                    break

            elif isoption(arg):
                optionname = arg#.split('-')[-1]
                plugin[optionname] = []
                plugin["_keyorder"].append(optionname)
                l = l[1:]
                while len(l):
                    arg = l[0]
                    if isplugin(arg) or isoption(arg):
                        break
                    plugin[optionname].append(evalifpossible(arg))
                    l = l[1:]

            else:
                plugin['main'].append(evalifpossible(arg))
                l = l[1:]
        return plugin, l

    # -----
    if len(sys.argv) == 1:
        return {}

    # -----
    D, remain = readplugin(sys.argv[1:])
    assert not len(remain)
    return D








if __name__ == "__main__":
    print (sys.argv[1:])
    D = readargv1()
    print (D.keys())
    #
    # print ("########MAIN")
    # print(D['main'])
    # print(D['_keyorder'])
    # print(D.keys())
    #
    # print("########option0")
    # print(D['option0'])
    #
    # print("########plugin1")
    # print(D['plugin1'])
    #
    # print("########plugin2")
    # print(D['plugin2'])
    # # print(D['plugin2'])
    #
