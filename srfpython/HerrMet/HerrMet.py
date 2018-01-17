#!/usr/bin/env python

import sys, matplotlib
if "-agg" in sys.argv[1:] or "--agg" in sys.argv[1:]:
    matplotlib.use('agg')

from tetedoeuf.utils.generic import readargv
from tetedoeuf.multipro.multipro8 import Job, MapAsync, MapSync

from srfpython.utils import depthspace, freqspace
from srfpython.inversion.metropolis2 import LogGaussND, LogUniND, metropolis
from srfpython.depthdisp.depthmodels import depthmodel, depthmodel1D, depthmodel_from_mod96
from srfpython.depthdisp.dispcurves import Claw
from srfpython.depthdisp.depthpdfs import dmstats1

import os, imp
import numpy as np
from displininv import Parameterizer, resultdisplay, Datacoder_log, makedatacoder, Theory, surf96reader
#from labex.ffmt.allfile import Allfile
#from labex.graphictools.gutils import  plt, gcf, gca, showme, grapher, values2colors, makecolorbar, chftsz, textonly
#from labex.inversion.metropolis2 import LogGaussND, LogUniND, metropolis
#from labex.processingtools.multipro7 import Job, MapAsync, MapSync
#from labex.dispersion.depthmodels import depthmodel, depthmodel1D, depthmodel_from_mod96, dmstats1, depthspace
#from labex.dispersion.displaws import Claw
#from labex.tools.stdout import waitbar
#from labex.signaltools.sutils import freqspace

#from labex.dispersion.VOLMAX.srfdis17 import groupbywtm, igroupbywtm, srfdis17
#from labex.dispersion.Herrmann0.dispersion import groupbywtm, igroupbywtm, dispersion as srfdis17

# -------------------------------------
def minmax(X):
    return min(X), max(X)
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
# -------------------------------------
# -------------------------------------
# -------------------------------------
def write_default_paramfile(nlayer, zbot, type = "mZVSPRRH", basedon=None, dvp=None, dvs=None, drh=None, dpr=None):
    """create a default parameter file to be customized by user"""
    # ------

    if np.all([_ is None for _ in dvs, dvp, drh, dpr]):
        which = None
    elif dvs is not None and dvp is None and drh is None and dpr is None:
        which = LogRhoM_DVS
    elif dvs is not None and dvp is not None and drh is not None and dpr is None:
        which = LogRhoM_DVPDVSDRH
    elif dvs is not None and dvp is not None and drh is not None and dpr is not None:
        which = LogRhoM_DVPDVSDRHDPR
    else:
        raise NotImplementedError('please specify either dvs alone, or dvp, dvs and drh, or dvp, dvs, drh and dpr')
    # ------
    def write_priortype_header(fid, dvp, dvs, drh):
        if which is None: pass
        elif which is LogRhoM_DVS:
            fid.write('#met PRIORTYPE = "DVS"\n')
            fid.write('#met DVSMIN = %f\n' % dvs[0])
            fid.write('#met DVSMAX = %f\n' % dvs[1])
        elif which is LogRhoM_DVPDVSDRH:
            fid.write('#met PRIORTYPE = "DVPDVSDRH"\n')
            fid.write('#met DVPMIN = %f\n' % dvp[0])
            fid.write('#met DVPMAX = %f\n' % dvp[1])
            fid.write('#met DVSMIN = %f\n' % dvs[0])
            fid.write('#met DVSMAX = %f\n' % dvs[1])
            fid.write('#met DRHMIN = %f\n' % drh[0])
            fid.write('#met DRHMAX = %f\n' % drh[1])
        elif which is LogRhoM_DVPDVSDRHDPR:
            fid.write('#met PRIORTYPE = "DVPDVSDRHDPR"\n')
            fid.write('#met DVPMIN = %f\n' % dvp[0])
            fid.write('#met DVPMAX = %f\n' % dvp[1])
            fid.write('#met DVSMIN = %f\n' % dvs[0])
            fid.write('#met DVSMAX = %f\n' % dvs[1])
            fid.write('#met DRHMIN = %f\n' % drh[0])
            fid.write('#met DRHMAX = %f\n' % drh[1])
            fid.write('#met DPRMIN = %f\n' % dpr[0])
            fid.write('#met DPRMAX = %f\n' % dpr[1])
        else:  raise Exception('programming error')

    # ------
    if type == "mZVSPRRH":
        if basedon is None:
            ztop = np.linspace(0, zbot, nlayer)
            ztopinf = -ztop[1:]  # deepest side
            ztopsup = -ztop[:-1] # shallow side
            ztopsup[0] -= 0.001
            vsinf = 0.1 * np.ones(nlayer)
            vssup = 3.5 * np.ones(nlayer)
            prinf = 1.6 * np.ones(nlayer) #r43 * np.ones(nlayer)
            prsup = 2.5 * np.ones(nlayer) #3.5 * np.ones(nlayer)
            rhinf = 2.2 * np.ones(nlayer)
            rhsup = 2.7 * np.ones(nlayer)
        else:
            b = depthmodel_from_mod96(basedon)
            #ztop = np.linspace(0., b.vp.z.max(), nlayer)
            ztop = depthspace(zbot, nlayer)
            b = b.interp(ztop, interpmethod="weightedstairs")
            ztopinf = -b.vp.z[1:]
            ztopsup = -b.vp.z[1:]
            vsinf = b.vs.values.copy()
            vssup = b.vs.values.copy()
            prinf = b.pr().values.copy()
            prsup = b.pr().values.copy()
            rhinf = b.rh.values.copy()
            rhsup = b.rh.values.copy()

        keys = ["-Z%d" % i for i in xrange(1, nlayer)] + \
               ["VS%d" % i for i in xrange(nlayer)] + \
               ["PR%d" % i for i in xrange(nlayer)] + \
               ["RH%d" % i for i in xrange(nlayer)]
        vinfs = np.concatenate((ztopinf, vsinf, prinf, rhinf))
        vsups = np.concatenate((ztopsup, vssup, prsup, rhsup))
        with open('_HerrMet.param', 'w') as fid:
            fid.write('#met TYPE = "mZVSPRRH"\n')
            fid.write('#met NLAYER = %d\n' % nlayer)
            write_priortype_header(fid, dvp, dvs, drh)
            fid.write('#fld KEY VINF VSUP\n')
            fid.write('#unt - - -\n')
            fid.write('#fmt %s  %f   %f\n')
            for k, vinf, vsup in zip(keys, vinfs, vsups):
                fid.write('%s %s %s\n' % (k, vinf, vsup))
        A = Allfile('_HerrMet.param')
        A.read()
        A.write() #to screen
        A.write('_HerrMet.param', fldorder = A._fldorder) #reformat file
    # ----------------------------
    elif type == "mZVSVPRH":
        if basedon is None:
            ztop = np.linspace(0, zbot, nlayer)
            ztopinf = -ztop[1:]  # deepest side
            ztopsup = -ztop[:-1] # shallow side
            ztopsup[0] -= 0.001
            vsinf = 0.1 * np.ones(nlayer)
            vssup = 3.5 * np.ones(nlayer)
            vpinf = 0.5 * np.ones(nlayer)
            vpsup = 6.5 * np.ones(nlayer)
            rhinf = 2.2 * np.ones(nlayer)
            rhsup = 2.7 * np.ones(nlayer)

        else:
            b = depthmodel_from_mod96(basedon)
            #ztop = np.linspace(0., b.vp.z.max(), nlayer)
            ztop = depthspace(zbot, nlayer)
            b = b.interp(ztop, interpmethod="weightedstairs")
            ztopinf = -b.vp.z[1:]
            ztopsup = -b.vp.z[1:]
            vsinf = b.vs.values.copy()
            vssup = b.vs.values.copy()
            vpinf = b.vp.values.copy()
            vpsup = b.vp.values.copy()
            rhinf = b.rh.values.copy()
            rhsup = b.rh.values.copy()


        keys = ["-Z%d" % i for i in xrange(1, nlayer)] + \
               ["VS%d" % i for i in xrange(nlayer)] + \
               ["VP%d" % i for i in xrange(nlayer)] + \
               ["RH%d" % i for i in xrange(nlayer)]

        vinfs = np.concatenate((ztopinf, vsinf, vpinf, rhinf))
        vsups = np.concatenate((ztopsup, vssup, vpsup, rhsup))

        with open('_HerrMet.param', 'w') as fid:
            fid.write('#met TYPE = "mZVSVPRH"\n')
            fid.write('#met NLAYER = %d\n' % nlayer)
            write_priortype_header(fid, dvp, dvs, drh)
            fid.write('#fld KEY VINF VSUP\n')
            fid.write('#unt - - -\n')
            fid.write('#fmt %s  %f   %f\n')
            for k, vinf, vsup in zip(keys, vinfs, vsups):
                fid.write('%s %s %s\n' % (k, vinf, vsup))
        A = Allfile('_HerrMet.param')
        A.read()
        A.write() #to screen
        A.write('_HerrMet.param', fldorder = A._fldorder) #reformat file
    # ----------------------------
    elif type == "mZVSPRzRHvp":
        if basedon is None:
            ztop = np.linspace(0, zbot, nlayer)
            ztopinf = -ztop[1:]  # deepest side
            ztopsup = -ztop[:-1] # shallow side
            ztopsup[0] -= 0.001
            vsinf = 0.1 * np.ones(nlayer)
            vssup = 3.5 * np.ones(nlayer)
        else:
            b = depthmodel_from_mod96(basedon)
            #ztop = np.linspace(0., b.vp.z.max(), nlayer)
            ztop = depthspace(zbot, nlayer)
            b = b.interp(ztop, interpmethod="weightedstairs")
            ztopinf = -b.vp.z[1:]
            ztopsup = -b.vp.z[1:]
            vsinf = b.vs.values.copy()
            vssup = b.vs.values.copy()

        keys = ["-Z%d" % i for i in xrange(1, nlayer)] + \
               ["VS%d" % i for i in xrange(nlayer)]

        vinfs = np.concatenate((ztopinf, vsinf))
        vsups = np.concatenate((ztopsup, vssup))

        with open('_HerrMet.param', 'w') as fid:
            fid.write('#met TYPE = "mZVSPRzRHvp"\n')
            fid.write('#met NLAYER = %d\n' % nlayer)
            fid.write('#met PRz  = "1.0335 * np.exp(-Z / 0.5408) + 1.7310  #some function of Z, Z is in km and growing downward"\n')
            fid.write('#met RHvp = "1.74 * VP ** 0.25 #some function of VP, VP is in km/s, RH is in g/cm3"\n')
            write_priortype_header(fid, dvp, dvs, drh)
            fid.write('#fld KEY VINF VSUP\n')
            fid.write('#unt - - -\n')
            fid.write('#fmt %s  %f   %f\n')
            for k, vinf, vsup in zip(keys, vinfs, vsups):
                fid.write('%s %s %s\n' % (k, vinf, vsup))
        A = Allfile('_HerrMet.param')
        A.read()
        A.write() #to screen
        A.write('_HerrMet.param', fldorder = A._fldorder) #reformat file
    # ----------------------------
    elif type == "mZVSPRzRHz":
        if basedon is None:
            ztop = np.linspace(0, zbot, nlayer)
            ztopinf = -ztop[1:]  # deepest side
            ztopsup = -ztop[:-1] # shallow side
            ztopsup[0] -= 0.001
            vsinf = 0.1 * np.ones(nlayer)
            vssup = 3.5 * np.ones(nlayer)
        else:
            b = depthmodel_from_mod96(basedon)
            #ztop = np.linspace(0., b.vp.z.max(), nlayer)
            ztop = depthspace(zbot, nlayer)
            b = b.interp(ztop, interpmethod="weightedstairs")
            ztopinf = -b.vp.z[1:]
            ztopsup = -b.vp.z[1:]
            vsinf = b.vs.values.copy()
            vssup = b.vs.values.copy()

        keys = ["-Z%d" % i for i in xrange(1, nlayer)] + \
               ["VS%d" % i for i in xrange(nlayer)]

        vinfs = np.concatenate((ztopinf, vsinf))
        vsups = np.concatenate((ztopsup, vssup))

        with open('_HerrMet.param', 'w') as fid:
            fid.write('#met TYPE = "mZVSPRzRHvp"\n')
            fid.write('#met NLAYER = %d\n' % nlayer)
            fid.write('#met PRz  = "1.0335 * np.exp(-Z / 0.5408) + 1.7310  #some function of Z, Z is in km and growing downward"\n')
            fid.write('#met RHz  = "Z * 0. + 2.67 #some function of Z, Z is in km and growing downward"\n')
            write_priortype_header(fid, dvp, dvs, drh)
            fid.write('#fld KEY VINF VSUP\n')
            fid.write('#unt - - -\n')
            fid.write('#fmt %s  %f   %f\n')
            for k, vinf, vsup in zip(keys, vinfs, vsups):
                fid.write('%s %s %s\n' % (k, vinf, vsup))
        A = Allfile('_HerrMet.param')
        A.read()
        A.write() #to screen
        A.write('_HerrMet.param', fldorder = A._fldorder) #reformat file
    else:
        raise NotImplementedError('no such parameter file type implemented %s' % type)
# -------------------------------------
# Parameterizers = object to convert an array of parameters into smth that can be understood by srfdis17
# -------------------------------------
class CommonParameterizer(Parameterizer):
    """default methods"""
    # ------------------
    def boundaries(self):
        """
        the lower model is not necessarily obtained when all parameters reach their lowest boundary
        """
        Ztopsup, VPlow, VSlow, RHlow = self.inv(self.MINF)
        Ztopinf, VPhgh, VShgh, RHhgh = self.inv(self.MSUP)
        z = np.sort(np.unique(np.concatenate((Ztopinf, Ztopsup))))

        #--------------------
        def f(Zinf, Zsup, V, which):
            v1 = depthmodel1D(Zinf, V).interp(z, interpmethod="stairs")
            v2 = depthmodel1D(Zsup, V).interp(z, interpmethod="stairs")
            return depthmodel1D(z, which(np.concatenate(([v1], [v2]), axis = 0), axis = 0))
        # --------------------
        vplow = f(Ztopinf, Ztopsup, VPlow, np.min)
        vphgh = f(Ztopinf, Ztopsup, VPhgh, np.max)

        vslow = f(Ztopinf, Ztopsup, VSlow, np.min)
        vshgh = f(Ztopinf, Ztopsup, VShgh, np.max)

        rhlow = f(Ztopinf, Ztopsup, RHlow, np.min)
        rhhgh = f(Ztopinf, Ztopsup, RHhgh, np.max)

        prlow = depthmodel1D(z, vphgh.values / vslow.values)
        prhgh = depthmodel1D(z, vplow.values / vshgh.values)

        return vplow, vphgh, vslow, vshgh, rhlow, rhhgh, prlow, prhgh

# -------------------------------------
class Parameterizer_mZVSPRRH(CommonParameterizer):
    """the parameterizer object links the model array (m) to a set of variables that can
       be understood by the theory function (ZTOP, VP, VS, RH)
    """
    def __init__(self, A):
        """Initiate from a Allfile object, assumes that read method has been called"""
        assert A.metadata['TYPE'] == "mZVSPRRH"
        assert np.all(A.data.VINF <= A.data.VSUP)
        if np.all(A.data.VINF == A.data.VSUP):
            print "Warning : all parameters are locked"
        # --------
        self.NLAYER = A.metadata['NLAYER']
        self.I = A.data.VINF < A.data.VSUP #index of dimension used, other parameters are kept constant
        # --------
        self.MDEFAULT = 0.5 * (A.data.VINF + A.data.VSUP) #used for dimensions that are not part of the model
        self.MMEAN = self.MDEFAULT[self.I]
        self.MINF  = A.data.VINF[self.I]
        self.MSUP  = A.data.VSUP[self.I]
        self.MSTD  = 0.5 * (A.data.VSUP - A.data.VINF)[self.I]
        # --------

    def boundaries(self):
        """the lower model is not necessarily obtained when all parameters reach their lowest boundary"""
        Ztopsup, VPlow, VSlow, RHlow = self.inv(self.MINF)
        Ztopinf, VPhgh, VShgh, RHhgh = self.inv(self.MSUP)
        print Ztopinf
        print Ztopsup
        z = np.sort(np.unique(np.concatenate((Ztopinf, Ztopsup))))
        PRlow = VPlow / VSlow #because it is the way it has been defined in self.inv
        PRhgh = VPhgh / VShgh
        del VPlow, VPhgh
        # --------------------
        def f(Zinf, Zsup, V, which):
            v1 = depthmodel1D(Zinf, V).interp(z, interpmethod="stairs")
            v2 = depthmodel1D(Zsup, V).interp(z, interpmethod="stairs")
            return depthmodel1D(z, which(np.concatenate(([v1], [v2]), axis=0), axis=0))

        # --------------------
        prlow = f(Ztopinf, Ztopsup, PRlow, np.min)
        prhgh = f(Ztopinf, Ztopsup, PRhgh, np.max)

        vslow = f(Ztopinf, Ztopsup, VSlow, np.min)
        vshgh = f(Ztopinf, Ztopsup, VShgh, np.max)

        rhlow = f(Ztopinf, Ztopsup, RHlow, np.min)
        rhhgh = f(Ztopinf, Ztopsup, RHhgh, np.max)

        vplow = depthmodel1D(z, prlow.values * vslow.values)
        vphgh = depthmodel1D(z, prhgh.values * vshgh.values)

        return vplow, vphgh, vslow, vshgh, rhlow, rhhgh, prlow, prhgh


            # ------------------
    def keys(self):
         k = ["-Z%d" % i for i in xrange(1, self.NLAYER)] + \
             ['VS%d' % i for i in xrange(self.NLAYER)] + \
             ['PR%d' % i for i in xrange(self.NLAYER)] + \
             ['RH%d' % i for i in xrange(self.NLAYER)]
         return np.asarray(k)[self.I]

    # ------------------
    def __call__(self, ZTOP, VP, VS, RH):
        """the method that is called to convert a real depth model into a model array
           called inside Theory.__call__ just before srfdis17
           please keep consistent with self.inv
        """
        m = np.concatenate((-1.0 * ZTOP[1:], VS, VP / VS, RH))[self.I]
        return m

    # ------------------
    def inv(self, m):
        """the method that is called to convert model array
           into a true depth model that can be understood by srfdis17
           called for plotting (i.e. convert a model back to something more physical)
           please keep consistent with self.__call__
        """
        M = self.MDEFAULT.copy()
        M[self.I] = m #overwrites default values with the one provided in m

        nlayer = self.NLAYER#(len(M) + 1) / 4
        ZTOP = np.concatenate(([0.], -1.0 * M[: 1 * nlayer - 1]))
        VS = M[1 * nlayer - 1: 2 * nlayer - 1]
        PR = M[2 * nlayer - 1: 3 * nlayer - 1]
        RH = M[3 * nlayer - 1: 4 * nlayer - 1]
        VP = PR * VS

        return ZTOP, VP, VS, RH
# -------------------------------------
class Parameterizer_mZVSVPRH(CommonParameterizer):
    """see Parameterizer_mZVSPRRH for doc"""
    def __init__(self, A):
        assert A.metadata['TYPE'] == "mZVSVPRH"
        assert np.all(A.data.VINF <= A.data.VSUP)
        if np.all(A.data.VINF == A.data.VSUP):
            print "Warning : all parameters are locked"
        # --------
        self.NLAYER = A.metadata['NLAYER']
        self.I = A.data.VINF < A.data.VSUP
        # --------
        self.MDEFAULT = 0.5 * (A.data.VINF + A.data.VSUP)      # value to use if not a parameter
        self.MMEAN = self.MDEFAULT[self.I]                     #
        self.MINF = A.data.VINF[self.I]                        # lower boundary for each parameter
        self.MSUP = A.data.VSUP[self.I]                        # upper boundary for each parameter
        self.MSTD = 0.5 * (A.data.VSUP - A.data.VINF)[self.I]  # markov proposal for each parameter
        # --------

    def boundaries(self):
        """
        the lower model is not necessarily obtained when all parameters reach their lowest boundary
        """
        Ztopsup, VPlow, VSlow, RHlow = self.inv(self.MINF)
        Ztopinf, VPhgh, VShgh, RHhgh = self.inv(self.MSUP)
        z = np.sort(np.unique(np.concatenate((Ztopinf, Ztopsup))))
        #--------------------
        def f(Zinf, Zsup, V, which):
            v1 = depthmodel1D(Zinf, V).interp(z, interpmethod="stairs")
            v2 = depthmodel1D(Zsup, V).interp(z, interpmethod="stairs")
            return depthmodel1D(z, which(np.concatenate(([v1], [v2]), axis = 0), axis = 0))
        # --------------------
        vplow = f(Ztopinf, Ztopsup, VPlow, np.min)
        vphgh = f(Ztopinf, Ztopsup, VPhgh, np.max)

        vslow = f(Ztopinf, Ztopsup, VSlow, np.min)
        vshgh = f(Ztopinf, Ztopsup, VShgh, np.max)

        rhlow = f(Ztopinf, Ztopsup, RHlow, np.min)
        rhhgh = f(Ztopinf, Ztopsup, RHhgh, np.max)

        prlow = depthmodel1D(z, vphgh.values / vslow.values)
        prhgh = depthmodel1D(z, vplow.values / vshgh.values)

        return vplow, vphgh, vslow, vshgh, rhlow, rhhgh, prlow, prhgh

    # ------------------
    def keys(self):
        k = ["-Z%d" % i for i in xrange(1, self.NLAYER)] + \
            ['VS%d' % i for i in xrange(self.NLAYER)] + \
            ['VP%d' % i for i in xrange(self.NLAYER)] + \
            ['RH%d' % i for i in xrange(self.NLAYER)]
        return np.asarray(k)[self.I]

    # ------------------
    def __call__(self, ZTOP, VP, VS, RH):
        m = np.concatenate((-1.0 * ZTOP[1:], VS, VP, RH))[self.I]
        return m

    # ------------------
    def inv(self, m):
        M = self.MDEFAULT.copy()
        M[self.I] = m  # overwrites default values with the one provided in m

        nlayer = self.NLAYER  # (len(M) + 1) / 4
        ZTOP = np.concatenate(([0.], -1.0 * M[: 1 * nlayer - 1]))
        VS = M[1 * nlayer - 1: 2 * nlayer - 1]
        VP = M[2 * nlayer - 1: 3 * nlayer - 1]
        RH = M[3 * nlayer - 1: 4 * nlayer - 1]

        return ZTOP, VP, VS, RH
# -------------------------------------
class Parameterizer_mZVSPRzRHvp(CommonParameterizer):
    """see Parameterizer_mZVSPRRH for doc"""

    def __init__(self, A):
        assert A.metadata['TYPE'] == "mZVSPRzRHvp"
        assert np.all(A.data.VINF <= A.data.VSUP)
        if np.all(A.data.VINF == A.data.VSUP):
            print "Warning : all parameters are locked"
        # --------
        self.NLAYER = A.metadata['NLAYER']
        self.I = A.data.VINF < A.data.VSUP
        # --------
        self.MDEFAULT = 0.5 * (A.data.VINF + A.data.VSUP)
        self.MMEAN = self.MDEFAULT[self.I]
        self.MINF = A.data.VINF[self.I]
        self.MSUP = A.data.VSUP[self.I]
        self.MSTD = 0.5 * (A.data.VSUP - A.data.VINF)[self.I]
        # --------
        # e.g. "1.0335 * np.exp(-Z / 0.5408) + 1.7310"
        self.PRzName=A.metadata['PRz'].split('#')[0].replace('np.', '')
        if "return" not in self.PRzName:
            self.PRz = string2func("import numpy as np\ndef PR(Z): return %s" % A.metadata['PRz'])
        else:
            self.PRz = string2func("import numpy as np\ndef PR(Z): %s" % A.metadata['PRz'])
            self.PRzName = ""
        # --------
        # e.g. "1.0335 * np.exp(-Z / 0.5408) + 1.7310"
        self.RHvpName = A.metadata['RHvp'].split('#')[0].replace('np.', '')
        self.RHvp = string2func("import numpy as np\ndef RH(VP): return %s" % A.metadata['RHvp'])
        # --------

    def boundaries(self):
        print "boundaries method not implemented for %s, using default one" % self.__class__.__name__
        return CommonParameterizer.boundaries(self)
    # ------------------
    def keys(self):
        k = ["-Z%d" % i for i in xrange(1, self.NLAYER)] + \
            ['VS%d' % i for i in xrange(self.NLAYER)]
        return np.asarray(k)[self.I]

    # ------------------
    def __call__(self, ZTOP, VP, VS, RH):
        """VP and RH are ignored since they are not parameters of the model"""
        m = np.concatenate((-1.0 * ZTOP[1:], VS))[self.I]
        return m

    # ------------------
    def inv(self, m):
        M = self.MDEFAULT.copy()
        M[self.I] = m  # overwrites default values with the one provided in m

        nlayer = self.NLAYER  # (len(M) + 1) / 4
        ZTOP = np.concatenate(([0.], -1.0 * M[: 1 * nlayer - 1]))
        VS = M[1 * nlayer - 1: 2 * nlayer - 1]

        #infer VP and RH from VS and input laws
        ZMID = np.concatenate((0.5 * (ZTOP[1:] + ZTOP[:-1]), [1.5 * ZTOP[-1]]))
        VP = VS * self.PRz(Z = ZMID)
        RH = self.RHvp(VP)

        return ZTOP, VP, VS, RH
# -------------------------------------
class Parameterizer_mZVSPRzRHz(CommonParameterizer):
    """see Parameterizer_mZVSPRRH for doc"""

    def __init__(self, A):
        assert A.metadata['TYPE'] == "mZVSPRzRHz"
        assert np.all(A.data.VINF <= A.data.VSUP)
        if np.all(A.data.VINF == A.data.VSUP):
            print "Warning : all parameters are locked"
        # --------
        self.NLAYER = A.metadata['NLAYER']
        self.I = A.data.VINF < A.data.VSUP
        # --------
        self.MDEFAULT = 0.5 * (A.data.VINF + A.data.VSUP)
        self.MMEAN = self.MDEFAULT[self.I]
        self.MINF = A.data.VINF[self.I]
        self.MSUP = A.data.VSUP[self.I]
        self.MSTD = 0.5 * (A.data.VSUP - A.data.VINF)[self.I]
        # --------
        # e.g. "1.0335 * np.exp(-Z / 0.5408) + 1.7310"
        self.PRzName=A.metadata['PRz'].split('#')[0].replace('np.', '')
        if "return" not in self.PRzName:
            self.PRz = string2func("import numpy as np\ndef PR(Z): return %s" % A.metadata['PRz'])
        else:
            self.PRz = string2func("import numpy as np\ndef PR(Z): %s" % A.metadata['PRz'])
            self.PRzName = ""
        # --------
        # e.g. "1.0335 * np.exp(-Z / 0.5408) + 1.7310"
        self.RHzName = A.metadata['RHz'].split('#')[0].replace('np.', '')
        if "return" not in self.RHzName:
            self.RHz = string2func("import numpy as np\ndef RH(Z): return %s" % A.metadata['RHz'])
        else:
            self.RHz = string2func("import numpy as np\ndef RH(Z): %s" % A.metadata['RHz'])
            self.RHzName = ""
        # --------

    def boundaries(self):
        print "boundaries method not implemented for %s, using default one" % self.__class__.__name__
        return CommonParameterizer.boundaries(self)
    # ------------------
    def keys(self):
        k = ["-Z%d" % i for i in xrange(1, self.NLAYER)] + \
            ['VS%d' % i for i in xrange(self.NLAYER)]
        return np.asarray(k)[self.I]

    # ------------------
    def __call__(self, ZTOP, VP, VS, RH):
        """VP and RH are ignored since they are not parameters of the model"""
        m = np.concatenate((-1.0 * ZTOP[1:], VS))[self.I]
        return m

    # ------------------
    def inv(self, m):
        M = self.MDEFAULT.copy()
        M[self.I] = m  # overwrites default values with the one provided in m

        nlayer = self.NLAYER  # (len(M) + 1) / 4
        ZTOP = np.concatenate(([0.], -1.0 * M[: 1 * nlayer - 1]))
        VS = M[1 * nlayer - 1: 2 * nlayer - 1]

        #infer VP and RH from VS and input laws
        ZMID = np.concatenate((0.5 * (ZTOP[1:] + ZTOP[:-1]), [1.5 * ZTOP[-1]]))
        VP = VS * self.PRz(Z = ZMID)
        RH = self.RHz(Z = ZMID)

        return ZTOP, VP, VS, RH
# -------------------------------------
# Prior probability density functions
# -------------------------------------
class DefaultLogRhoM(LogUniND):
    """build the prior pdf, include constraints on the model paramters and eventually on some relations between them
    """
    k, nanbehavior = 1000., 2
    # -------
    def __init__(self, parameterizer):
        """no more constraints than parameters"""
        LogUniND.__init__(self,
                      vinfs=parameterizer.MINF,
                      vsups=parameterizer.MSUP,
                      k=self.k,
                      nanbehavior=self.nanbehavior)
# -------------------------------------
class LogRhoM_DVS(DefaultLogRhoM):
    """add new constraitns on the vs offsets on interfaces"""
    # -------
    def __init__(self, parameterizer, dvsmin, dvsmax):
        assert dvsmin < dvsmax
        self.p = parameterizer
        NLAYER = self.p.NLAYER
        MINF = np.concatenate((self.p.MINF, np.ones(NLAYER - 1) * dvsmin))
        MSUP = np.concatenate((self.p.MSUP, np.ones(NLAYER - 1) * dvsmax))
        LogUniND.__init__(self,
                      vinfs=MINF,
                      vsups=MSUP,
                      k=self.k,
                      nanbehavior=self.nanbehavior)
    # -------
    def __call__(self, m):
        _, _, VS, _ = self.p.inv(m)
        DVS = VS[1:] - VS[:-1]
        extended_m = np.concatenate((m, DVS))
        return LogUniND.__call__(self, extended_m)
# -------------------------------------
class LogRhoM_DVPDVSDRH(DefaultLogRhoM):
    """add new constraitns on the vs offsets on interfaces"""
    # -------
    def __init__(self, parameterizer, dvpmin, dvpmax, dvsmin, dvsmax, drhmin, drhmax):
        assert dvsmin < dvsmax
        assert dvpmin < dvpmax
        assert drhmin < drhmax
        self.p = parameterizer
        NLAYER = self.p.NLAYER
        MINF = np.concatenate((self.p.MINF,
                               np.ones(NLAYER - 1) * dvpmin,
                               np.ones(NLAYER - 1) * dvsmin,
                               np.ones(NLAYER - 1) * drhmin))
        MSUP = np.concatenate((self.p.MSUP,
                               np.ones(NLAYER - 1) * dvpmax,
                               np.ones(NLAYER - 1) * dvsmax,
                               np.ones(NLAYER - 1) * drhmax))
        LogUniND.__init__(self,
                      vinfs=MINF,
                      vsups=MSUP,
                      k=self.k,
                      nanbehavior=self.nanbehavior)
    # -------
    def __call__(self, m):
        _, VP, VS, RH = self.p.inv(m)
        DVP = VS[1:] - VS[:-1]
        DVS = VP[1:] - VP[:-1]
        DRH = RH[1:] - RH[:-1]
        extended_m = np.concatenate((m, DVP, DVS, DRH))
        return LogUniND.__call__(self, extended_m)
# -------------------------------------
class LogRhoM_DVPDVSDRHDPR(DefaultLogRhoM):
    """add new constraitns on the vs offsets on interfaces"""

    # -------
    def __init__(self, parameterizer, dvpmin, dvpmax, dvsmin, dvsmax, drhmin, drhmax, dprmin, dprmax):
        assert dvsmin < dvsmax
        assert dvpmin < dvpmax
        assert drhmin < drhmax
        assert dprmin < dprmax
        self.p = parameterizer
        NLAYER = self.p.NLAYER
        MINF = np.concatenate((self.p.MINF,
                               np.ones(NLAYER - 1) * dvpmin,
                               np.ones(NLAYER - 1) * dvsmin,
                               np.ones(NLAYER - 1) * drhmin,
                               np.ones(NLAYER - 1) * dprmin))
        MSUP = np.concatenate((self.p.MSUP,
                               np.ones(NLAYER - 1) * dvpmax,
                               np.ones(NLAYER - 1) * dvsmax,
                               np.ones(NLAYER - 1) * drhmax,
                               np.ones(NLAYER - 1) * dprmax))
        LogUniND.__init__(self,
                          vinfs=MINF,
                          vsups=MSUP,
                          k=self.k,
                          nanbehavior=self.nanbehavior)

    # -------
    def __call__(self, m):
        _, VP, VS, RH = self.p.inv(m)
        DVP = VS[1:] - VS[:-1]
        DVS = VP[1:] - VP[:-1]
        DRH = RH[1:] - RH[:-1]
        PR = VP / VS
        DPR = PR[1:] - PR[:-1]
        extended_m = np.concatenate((m, DVP, DVS, DRH, DPR))
        return LogUniND.__call__(self, extended_m)

# -------------------------------------
# -------------------------------------
# -------------------------------------
def load_paramfile(paramfile):
    """initiate one of the parameterizer and prior pdf according to the param file"""
    A = Allfile(paramfile)
    A.read()

    # ------------------------
    if A.metadata['TYPE'] == "mZVSVPRH":      p = Parameterizer_mZVSVPRH(A)
    elif A.metadata['TYPE'] == "mZVSPRRH":    p = Parameterizer_mZVSPRRH(A)
    elif A.metadata['TYPE'] == "mZVSPRzRHvp": p = Parameterizer_mZVSPRzRHvp(A)
    elif A.metadata['TYPE'] == "mZVSPRzRHz":  p = Parameterizer_mZVSPRzRHz(A)
    else: raise Exception('could not load %s (TYPE not recognized)' % paramfile)

    # ------------------------
    if not "PRIORTYPE" in A.metadata.keys():
        logRHOM = DefaultLogRhoM(p)
    elif A.metadata['PRIORTYPE'] == "DVS":
        logRHOM = LogRhoM_DVS(p,
                    dvsmin=A.metadata['DVSMIN'],
                    dvsmax=A.metadata['DVSMAX'])
    elif A.metadata['PRIORTYPE'] == "DVPDVSDRH":
        logRHOM = LogRhoM_DVPDVSDRH(p,
                    dvpmin=A.metadata['DVPMIN'],
                    dvpmax=A.metadata['DVPMAX'],
                    dvsmin=A.metadata['DVSMIN'],
                    dvsmax=A.metadata['DVSMAX'],
                    drhmin=A.metadata['DRHMIN'],
                    drhmax=A.metadata['DRHMAX'])
    elif A.metadata['PRIORTYPE'] == "DVPDVSDRHDPR":
        logRHOM = LogRhoM_DVPDVSDRHDPR(p,
                    dvpmin=A.metadata['DVPMIN'],
                    dvpmax=A.metadata['DVPMAX'],
                    dvsmin=A.metadata['DVSMIN'],
                    dvsmax=A.metadata['DVSMAX'],
                    drhmin=A.metadata['DRHMIN'],
                    drhmax=A.metadata['DRHMAX'],
                    dprmin=A.metadata['DPRMIN'],
                    dprmax=A.metadata['DPRMAX'])
    else:
        raise Exception('could not load %s (PRIORTYPE not recognized)' % paramfile)


    # ------------------------
    print "parameter type : ", p.__class__.__name__
    print "prior type     : ", logRHOM.__class__.__name__
    return p, logRHOM
# -------------------------------------
def readHerrmetout_serial(f):
    with open(f, 'r') as fid:
        nfield = None
        while True:
            l = fid.readline().strip()
            if l == "": break
            if l == "\n" or l.startswith("#"): continue
            l = l.strip('\n').split()
            if nfield is None: nfield = len(l)
            elif not len(l) == nfield:
                print l
                break #line is not full, a run is probably in progress
            chainid, weight, nlayer = np.asarray(l[:3], int)
            llk = float(l[3])
            lmod = np.concatenate(([0.], np.asarray(l[4:4 + 4 * nlayer - 1], float)))

            ldat = l[4 + 4 * nlayer - 1:]
            ndat = len(ldat) / 5

            ztop = lmod[:nlayer]
            vp = lmod[nlayer: 2 * nlayer]
            vs = lmod[2 * nlayer: 3 * nlayer]
            rh = lmod[3 * nlayer: 4 * nlayer]
            waves = ldat[:ndat]
            types = ldat[ndat:2 * ndat]
            modes = np.array(ldat[2 * ndat:3 * ndat], int)
            freqs = np.array(ldat[3 * ndat:4 * ndat], float)
            values = np.array(ldat[4 * ndat:5 * ndat], float)
            yield chainid, weight, llk, (ztop, vp, vs, rh), (waves, types, modes, freqs, values)
# -------------------------------------
def readHerrmetout(f, **mapkwargs):
    def gen():
        with open(f, 'r') as fid:
            nfield = None
            while True:
                l = fid.readline().strip()
                if l == "": break
                if l == "\n" or l.startswith("#"): continue
                l = l.strip('\n').split()
                if nfield is None:
                    nfield = len(l)
                elif not len(l) == nfield:
                    print l
                    break  # line is not full, a run is probably in progress
                yield Job(l)
    def fun(l):
        chainid, weight, nlayer = np.asarray(l[:3], int)
        llk = float(l[3])
        lmod = np.concatenate(([0.], np.asarray(l[4:4 + 4 * nlayer - 1], float)))

        ldat = l[4 + 4 * nlayer - 1:]
        ndat = len(ldat) / 5

        ztop = lmod[:nlayer]
        vp = lmod[nlayer: 2 * nlayer]
        vs = lmod[2 * nlayer: 3 * nlayer]
        rh = lmod[3 * nlayer: 4 * nlayer]
        waves = ldat[:ndat]
        types = ldat[ndat:2 * ndat]
        modes = np.array(ldat[2 * ndat:3 * ndat], int)
        freqs = np.array(ldat[3 * ndat:4 * ndat], float)
        values = np.array(ldat[4 * ndat:5 * ndat], float)
        return  chainid, weight, llk, (ztop, vp, vs, rh), (waves, types, modes, freqs, values)
    with MapAsync(fun, gen(), **mapkwargs) as ma:
        for _, ans, _, _ in ma:
            #print (j[1] - j[0]) / (g[1] - g[0])
            yield ans
# -------------------------------------
def readHerrmetout_1(f, top=None, topstep=1, **mapkwargs):
    chainids, weights, llks, ms, ds = zip(*list(readHerrmetout(f, **mapkwargs)))
    chainids, weights, llks = [np.asarray(_) for _ in chainids, weights, llks]
    if top is not None:
        I = np.argsort(llks)[::-1][:top][::topstep]
    else: #means all data
        I = np.argsort(llks)[::-1]
    return chainids[I], weights[I], llks[I], [ms[i] for i in I], [ds[i] for i in I]
# -------------------------------------
def readHerrmetout_2(f, **mapkwargs):
    with open(f, 'r') as fid:
        nmodels = 1
        while True:
            l = fid.readline()
            if l == "": break
            if l.startswith('#'): continue
            nmodels += 1
    # -------------------
    g = readHerrmetout(f, **mapkwargs)
    chainid, weight, llk, (ztop, vp, vs, rh), (waves, types, modes, freqs, values) = g.next()
    ZTOP = np.zeros((nmodels, len(ztop)), float) * np.nan
    VP = np.zeros((nmodels, len(vp)), float) * np.nan
    VS   = np.zeros((nmodels, len(vs)), float) * np.nan
    RH = np.zeros((nmodels, len(rh)), float) * np.nan
    WEIGHTS = np.zeros(nmodels, int)
    LLKS = np.zeros(nmodels, int) * np.nan
    # -------------------
    ZTOP[0, :], VP[0, :], VS[0, :], RH[0,:], WEIGHTS[0], LLKS[0] = ztop, vp, vs, rh, weight, llk
    # -------------------
    for n, (chainid, weight, llk, (ztop, vp, vs, rh), _) in enumerate(readHerrmetout(f)):
        ZTOP[n+1, :], \
        VP[n+1, :], \
        VS[n+1, :], \
        RH[n+1, :], \
        WEIGHTS[n+1],\
        LLKS = ztop, vp, vs, rh, weight, llk
    # -------------------
    return ZTOP, VP, VS, RH, WEIGHTS, LLKS
# -------------------------------------
def overdisp(ms, overwaves, overtypes, overmodes, overfreqs, **mapkwargs):
    """extrapolate dispersion curves"""
    def fun(mms):
        ztop, vp, vs, rh = mms
        try:
            overvalues = srfdis17(ztop, vp, vs, rh, overwaves, overtypes, overmodes, overfreqs, h = 0.005, dcl = 0.005, dcr = 0.005)
        except KeyboardInterrupt: raise
        except Exception as e:
            h = ztop[1:] - ztop[:-1]
            # assume failuer was caused by rounding issues
            h[h <= 0.001] = 0.001001
            ztop = np.concatenate(([0.], h.cumsum()))
            try: #again
                overvalues = srfdis17(ztop, vp, vs, rh, overwaves, overtypes, overmodes, overfreqs, h=0.005, dcl=0.005,
                                      dcr=0.005)
            except KeyboardInterrupt: raise
            except Exception as giveup:
                overvalues = np.nan * np.ones(len(overwaves))
        return mms, overvalues

    with MapSync(fun, (Job(mms) for mms in ms), **mapkwargs) as ma:
        wb = waitbar('overdisp')
        Njobs = len(ms) - 1.
        for jobid, (mms, overvalues), _, _ in ma:
            wb.refresh(jobid / Njobs)
            dds = (overwaves, overtypes, overmodes, overfreqs, overvalues)
            yield mms, dds
        wb.close()
    print
# -------------------------------------


# -------------------------------------
version = "4.1"
default_mode = "append"
default_nchain = 12
default_nkeep = 100
default_top = 100
default_topstep = 1
default_parameterization_list = ['mZVSPRRH', 'mZVSVPRH', 'mZVSPRzRHvp', 'mZVSPRzRHz']
default_parameterization = default_parameterization_list[0]
mapkwargs = {} #keyword arguments for every parallelized process
# -------------------------------------
autorizedkeys = \
    ["w", "taskset", "agg", "lowprio",
     "help", "h",
     "example", "ex",
     "param", "basedon", "t", "dvp", "dvs", "drh", "growing", "op",
     "target", "resamp", "lunc", "unc", "ot",
     "run", "nchain", "nkeep",
     "extract",
     "disp", "best", "overdisp", "range", "png", "m96",
     "test", "sltz", "ritt"]
# -------------------------------------
help = '''HerrMet V{version}
# --------------------------
-w                     change number of virtual workers for all parallelized processes
-taskset               change job affinity for all parallelized processes, e.g. "0-4"
-agg                   use agg backend (no display)
-lowprio               run processes with low priority
# --------------------------
--help, -h             display this help message, and quit
--example, -e          display an example of script, and quit
--param      int float generate a template parameter file to custom, need the number of layers and bottom depth in km
    -basedon filename  build parametrization based on an existing mod96 file, if not specified, I take fixed values
    -t       typename  parameterization type to use ({default_parameterization_list}), default {default_parameterization}
                       mZVSPRRH = parameterize with depth interface, 
                                  VS in each layer, VP/VS in each layer, Density in each layer
                       mZVSVPRH = parameterize with depth interface,   
                                  VP in each layer, VP/VS in each layer, Density in each layer
                       mZVSPRzRHvp = parameterize with depth interface, 
                                  use fixed relations for VP/VS versus depth and Density versus VP
                      mZVSPRzRHz = parameterize with depth interface, 
                                  use fixed relations for VP/VS versus depth and Density versus depth 
    -dvp     f f       add prior constraint on the vp offset between layers, provide minimum and maximum value, km/s
    -dvs     f f       add prior constraint on the vs offset between layers, provide minimum and maximum value, km/s
    -drh     f f       add prior constraint on the density offset between layers, provide minimum and maximum value, g/cm3
    -dpr     f f       add prior constraint on the vp/vs offset between layers, provide minimum and maximum value
    -growing           shortcut for -dvp 0. 5. -dvs 0. 5. -drh 0. 5. -dpr -5. 0.
    -op                overwrite _HerrMet.param if exists
--target     filename  set the target dispersion curve from surf96
    -resamp  f f i s   resample the dispersion curve, 
                       needs fmin(Hz), fmax(Hz), nfreq, fscale(flin, plin or log)
    -lunc    float     set constant uncertainty in log domain (value x lunc)
    -unc     float     set constant uncertainty in linear domain 
    -ot                overwrite _HerrMet.target if exists
--run        mode      start inversion, mode is append or restart
    -nchain  int       number of chains to use, default {default_nchain}
    -nkeep   int       number of models to keep per chain, default {default_nkeep}
    -w                 see above, controls the max number of chains to be run simultaneously
--extract   [i] [i]    extract posterior distribution on the models, save them as mod96files
                       first argument = number of best models to use/display, default {default_top}
                       second argument = step between them, default {default_topstep}
--disp [i] [i]         display param, target, and run outputs if exists
                       first argument = number of best models to use/display, default {default_top}
                       second argument = step between them, default {default_topstep}
    -best/-overdisp    show the best models on the figure, use overdisp instead of best to recompute dispersion curves with higher resolution
    -range             compute and show the statistics for the selected models, use --extract for saving
    -png               save figure as pngfile instead of displaying it on screen
    -m96 file(s)       append depth model(s) to the plot from mod96 file(s)    
--test                 testing option
'''.format(
    version=version,
    default_top=default_top,
    default_nchain=default_nchain,
    default_nkeep=default_nkeep,
    default_topstep=default_topstep,
    default_parameterization_list=default_parameterization_list,
    default_parameterization=default_parameterization
    )
# -------------------------------------
example="""
# get target, resample it between 0.2-1.5 Hz with 15 samples spaced logarithmically in period domain
# adjust uncertainties to 0.1 in logaritmic domain, overwrite target if exists (_HerrMet.target) 
# and display it
HerrMet --target /path/to/my/data/file.surf96 \\
            -resamp 0.2 1.5 15 plog \\
            -lunc 0.1 \\
            -ot \\
            --disp

# build parameter file from existing depthmodel, use 7 layers, use parametrization mZVSPRRH, 
# require vp, vs and density to be growing
# overwrite paramfile if exists (_HerrMet.param) and display
HerrMet --param 7 \\
            -basedon /path/to/my/depthmodel.mod96 \\
            -t  mZVSPRRH \\
            -growing \\
            -op \\
            --disp

# >> now edit _HerrMet.param and custom it, check with HerrMet --disp

# run inversion with 12 chains, keep 1000 models each, run on 24 virtual threads
HerrMet -w 24 \\
        --run restart \\
            -nchain 12 -nkeep 1000 
        
# display best 1000 models, recompute the best disp curves with higher resolution
# compute median and percentiles over these 1000 models
# save as png file, use a non display backend 
HerrMet -agg \\
        --disp 1000 \\
            -overdisp \\
            -range \\
            -png
"""


# -------------------------------------
if __name__ == "__main__":
    if len(sys.argv) == 1:
        print help
        sys.exit()
    argv = readargv()
    # nworkers = int(argv['w'][0]) if "w" in argv.keys() else default_nworkers
    # -------------------------------------
    # prevent typos in arguments, keep the autorizedkeys list up to date
    for k in argv.keys():
        if k not in autorizedkeys:
            if k.startswith('_'): continue # private keys
            raise Exception('keyword %s is not recognized' % k)
    # -------------------------------------
    if "w" in argv.keys():
        mapkwargs["Nworkers"] = int(argv['w'][0])
    if "taskset" in argv.keys():
        mapkwargs["Taskset"] = argv['taskset'][0]
        os.system('taskset -pc %s %d' % (argv['taskset'][0], os.getpid()))
    if "lowprio" in argv.keys():
        mapkwargs["LowPriority"] = True
    # -------------------------------------
    if "h" in argv.keys() or "help" in argv.keys():
        print help
        sys.exit()
    # -------------------------------------
    if "e" in argv.keys() or "example" in argv.keys():
        print example
        sys.exit()
    # -------------------------------------
    elif "v" in argv.keys() or "version" in argv.keys():
        print "version : %s" % version
        sys.exit()
    # -------------------------------------
    if "clean" in argv.keys():
        os.system('rm -f ./_HerrMet.param ./_HerrMet.target ./_HerrMet.run ./_HerrMet.p*.mod96')
        # sys.exit()
    # -------------------------------------
    if "param" in argv.keys():
        if "op" not in argv.keys():
            assert not os.path.exists('_HerrMet.param')

        nlayer = int(argv["param"][0])
        zbot   = float(argv["param"][1])
        type = argv['t'][0] if "t" in argv.keys() else default_parameterization
        basedon = argv['basedon'][0] if "basedon" in argv.keys() else None
        if "growing" in argv.keys():
            assert "dvs" not in argv.keys() #not compatible with -growing
            assert "dvp" not in argv.keys() #not compatible with -growing
            assert "drh" not in argv.keys() #not compatible with -growing
            assert "dpr" not in argv.keys() # not compatible with -growing
            dvp = 0., 5.
            dvs = 0., 5.
            drh = 0., 5.
            dpr = -5., 0.
        else:
            dvp = minmax(np.asarray(argv['dvp'], float)) if "dvp" in argv.keys() else None
            dvs = minmax(np.asarray(argv['dvs'], float)) if "dvs" in argv.keys() else None
            drh = minmax(np.asarray(argv['drh'], float)) if "drh" in argv.keys() else None
            dpr = minmax(np.asarray(argv['dpr'], float)) if "dpr" in argv.keys() else None


        if not type in default_parameterization_list:
            raise Exception('please pick one type in %s' % str(default_parameterization_list))
        write_default_paramfile(nlayer, zbot, type=type, basedon=basedon, dvp=dvp, dvs=dvs, drh=drh, dpr=dpr)
        print "please customize _HerrMet.param, do not change line orders and metadata"
        print "use option --disp to see the depth boudaries"
        # sys.exit()
    # -------------------------------------
    if "target" in argv.keys():
        if "ot" not in argv.keys():
            assert not os.path.exists('_HerrMet.target')

        s96 = argv["target"][0]
        s = surf96reader(s96)
        # -------------------
        if "resamp" in argv.keys():
            news = s.copy()
            news.accept(np.zeros(len(news), bool)) #clear all entries
            newf = freqspace(freqmin=float(argv["resamp"][0]),
                             freqmax=float(argv["resamp"][1]),
                             nfreq=int(argv["resamp"][2]),
                             scale=argv["resamp"][3])
            for law in s.get_all():
                law.set_extrapolationmode(1)
                stdlaw = Claw(freq = law.freq, value = law.dvalue, extrapolationmode = 0)

                newvalues = law(newf)
                newdvalues = stdlaw(newf)

                I = ~np.isnan(newvalues)
                if I.any():
                    N = I.sum()
                    news.data['WAVE'] = np.concatenate((news.data['WAVE'], np.array([law.wave]).repeat(N)))
                    news.data['TYPE'] = np.concatenate((news.data['TYPE'], np.array([law.type]).repeat(N)))
                    news.data['MODE'] = np.concatenate((news.data['MODE'], np.array([law.mode]).repeat(N)))
                    news.data['PERIOD'] = np.concatenate((news.data['PERIOD'], 1. / newf[I]))
                    news.data['VALUE'] = np.concatenate((news.data['VALUE'], newvalues[I]))
                    news.data['DVALUE'] = np.concatenate((news.data['DVALUE'], newdvalues[I]))
            s = news
            print news
        # -------------------
        if "lunc" in argv.keys():
            # set uncertainties to constant in log domain
            lunc = float(argv["lunc"][0])
            s.data['DVALUE'] = s.data['VALUE'] * lunc
            s.write96('_HerrMet.target')
        elif "unc" in argv.keys():
            # set uncertainties to constant in lin domain
            unc = float(argv["unc"][0])
            s.data['DVALUE'] = unc
        # -------------------
        s.write96('_HerrMet.target')
        print "please only datapoints to invert in _HerrMet.target"
        print "use option --disp to see the target data"
        # sys.exit()
    # -------------------------------------
    if "run" in argv.keys():
        mode = argv['run'][0]
        assert mode in ['append', 'restart']

        Nchain = int(argv['nchain'][0]) if "nchain" in argv.keys() else default_nchain
        Nkeep  = int(argv['nkeep'][0]) if "nkeep" in argv.keys() else default_nkeep

        # ------
        p, logRHOM = load_paramfile('_HerrMet.param')
        # ------
        d = makedatacoder("_HerrMet.target", which=Datacoder_log)  # datacoder based on observations
        dobs, CDinv = d.target()
        duncs = CDinv ** -.5
        ND = len(dobs)
        dinfs = d(0.1 * np.ones_like(d.values))
        dsups = d(3.5 * np.ones_like(d.values))
        logRHOD = LogGaussND(dobs, duncs, dinfs, dsups, k=1000., nanbehavior=1)
        # ------
        G = Theory(parameterizer=p, datacoder=d)
        # ------
        # rd = resultdisplay1(targetfile="_HerrMet.target")
        # rd.plotmodel(alpha = 1.0, color = "r", linewidth = 3, *p.inv(p.MINF))
        # rd.plotmodel(alpha = 1.0, color="r", linewidth=3, *p.inv(p.MSUP))
        # rd.plotdisp(d.waves, d.types, d.modes, d.freqs, d.inv(dobs), dvalues = d.dvalues, alpha = 1.0, color="k", linewidth=3)
        # showme(False)
        # ---------------------------------
        def gen():
            for nchain in xrange(Nchain):
                M0 = np.random.rand(len(p.MINF)) * (p.MSUP - p.MINF) + p.MINF
                MSTD = p.MSTD
                yield Job(nchain, M0, MSTD, nkeep=Nkeep)


        def fun(worker, chainid, M0, MSTD, nkeep):
            models, datas, weights, llks = metropolis(M0, MSTD, G, ND, logRHOD, logRHOM,
                  nkeep=nkeep,
                  normallaw=worker.randn,
                  unilaw=worker.rand,
                  chainid=chainid,
                  HL=10,
                  IK0=0.25,
                  MPMIN=1.e-6,
                  MPMAX=1e6,
                  adjustspeed=0.3,
                  nofail=True,
                  debug=False)

            I = np.any(~np.isnan(datas), axis=1)
            models, datas, weights, llks = models[I, :], datas[I, :], weights[I], llks[I]

            return chainid, models, datas, weights, llks
        # ---------------------------------
        with MapAsync(fun, gen(), **mapkwargs) as ma, open('_HerrMet.run', 'w' if mode == "restart" else "a") as fid:
            fid.write('#CHAINID WEIGHT NLAYER LLK ZTOP[1:] VP VS RH WAVES TYPES MODES FREQS VALUES\n')
            for jobid, (chainid, models, datas, weights, llks), _, _ in ma:
                for mdl, dat, wgt, llk in zip(models, datas, weights, llks):
                    ztop, vp, vs, rh = p.inv(mdl)
                    values = d.inv(dat)
                    nlayer = len(ztop)
                    fid.write("%d %d %d %f %s %s %s %s %s %s %s %s %s\n" %
                              (chainid, wgt, nlayer, llk,
                               tostr(ztop[1:], "%.4f"),
                               tostr(vp, "%.3f"),
                               tostr(vs, "%.3f"),
                               tostr(rh, "%.3f"),
                               tostr(d.waves, "%s"),
                               tostr(d.types, "%s"),
                               tostr(d.modes, "%d"),
                               tostr(d.freqs, "%.4f"),
                               tostr(values, "%4f")))

        # sys.exit()

    # -------------------------------------
    if "extract" in  argv.keys():

        top = int(argv['extract'][0]) if argv['extract'] is not None else default_top
        topstep = int(argv['extract'][1]) if argv['extract'] is not None and len(argv['extract']) >= 2 else default_topstep

        chainids, weights, llks, ms, ds = readHerrmetout_1('_HerrMet.run', top=top, topstep=topstep)
        dms, wgts = [], []
        for weight, (ztop, vp, vs, rh) in zip(weights, ms):  # readHerrmetout("_HerrMet.run"):
            dm = depthmodel(depthmodel1D(ztop, vp),
                            depthmodel1D(ztop, vs),
                            depthmodel1D(ztop, rh))
            dms.append(dm)
            wgts.append(weight)
        for p, (vppc, vspc, rhpc, prpc) in dmstats1(dms,
                percentiles=[0.16, 0.5, 0.84],
                Ndepth=100,
                Nvalue=100,
                weights=wgts, **mapkwargs):
            try:
                dmout = depthmodel(vppc, vspc, rhpc)
                dmout.write96('_HerrMet.p%.2f.mod96' % p, overwrite=True)
            except KeyboardInterrupt: raise
            except Exception as e:
                print "Error", str(e)
    # -------------------------------------
    if "disp" in argv.keys():
        assert not ("best" in argv.keys() and "overdisp" in argv.keys()) #options are not compatible

        top = int(argv['disp'][0]) if argv['disp'] is not None else default_top
        topstep = int(argv['disp'][1]) if argv['disp'] is not None and len(argv['disp']) >= 2 else default_topstep

        # ------ Display the target data if exists
        if os.path.exists("_HerrMet.target"):
            rd = resultdisplay1(targetfile="_HerrMet.target")
            d = makedatacoder("./_HerrMet.target", which=Datacoder_log)  # datacoder based on observations
            dobs, _ = d.target()
        else:
            rd = resultdisplay1()
            print "call option --target to see the target dispersion curves"

        # ------ Display run results if exist
        if ("best" in argv.keys() or "range" in argv.keys() or "overdisp" in argv.keys()) and os.path.exists('_HerrMet.run'):

            chainids, weights, llks, ms, ds = readHerrmetout_1('_HerrMet.run', top=top, topstep=topstep)
            vmin, vmax, cmap = llks.min(), llks.max(),   plt.cm.gray  #plt.cm.jet# plt.cm.gray #
            colors = values2colors(llks, vmin=vmin, vmax=vmax, cmap=cmap)


            if "best" in argv.keys():
                for i in range(len(llks))[::-1]:
                    ztop, vp, vs, rh = ms[i]
                    dm = depthmodel(depthmodel1D(ztop, vp),
                                    depthmodel1D(ztop, vs),
                                    depthmodel1D(ztop, rh))
                    # dm.write96('M%010.0f.mod' % i, overwrite = True)
                    rd.plotmodel(color=colors[i], alpha=1.0, linewidth=3, *ms[i])
                    rd.plotdisp(color=colors[i], alpha=1.0, linewidth=3, *ds[i])

                cb = makecolorbar(vmin=vmin, vmax=vmax, cmap=cmap)
                cax = rd.fig.add_axes((0.78, 0.1, 0.005, 0.3))
                rd.fig.colorbar(cb, cax=cax, label="log likelyhood")

            elif "overdisp" in argv.keys():
                """note : recomputing dispersion with another frequency array might
                          result in a completely different dispersion curve in case
                          of root search failure """
                waves, types, modes, freqs, _ = ds[0]
                overwaves, overtypes, overmodes, _, _ = zip(*list(groupbywtm(waves, types, modes, freqs, np.arange(len(freqs)), None, True)))
                overfreqs = [freqspace(0.6 * min(freqs), 1.4 * max(freqs), 100, "plog") for _ in xrange(len(overwaves))]
                overwaves, overtypes, overmodes, overfreqs = igroupbywtm(overwaves, overtypes, overmodes, overfreqs)
                for clr, (mms, dds) in zip(colors[::-1], overdisp(ms[::-1], overwaves, overtypes, overmodes, overfreqs, **mapkwargs)):
                    rd.plotmodel(color=clr, alpha=1.0, linewidth=3, *mms)
                    try:
                        rd.plotdisp(color=clr, alpha=1.0, linewidth=3, *dds)
                    except KeyboardInterrupt: raise
                    except Exception as e:
                        print "Error : could not plot dispersion curve (%s)" % str(e)

                cb = makecolorbar(vmin=vmin, vmax=vmax, cmap=cmap)
                cax = rd.fig.add_axes((0.78, 0.1, 0.005, 0.3))
                rd.fig.colorbar(cb, cax=cax, label="log likelyhood")

            if "range" in argv.keys():
                dms, wgts = [], []
                for weight, (ztop, vp, vs, rh) in zip(weights, ms):#readHerrmetout("_HerrMet.run"):
                    dm = depthmodel(depthmodel1D(ztop, vp),
                                    depthmodel1D(ztop, vs),
                                    depthmodel1D(ztop, rh))
                    dms.append(dm)
                    wgts.append(weight)
                for p, (vppc, vspc, rhpc, prpc) in dmstats1(dms,
                            percentiles=[0.16, 0.5, 0.84],
                            Ndepth=100,
                            Nvalue=100,
                            weights=wgts, **mapkwargs):
                    # for _ in vppc, vspc, rhpc, prpc:
                    #     _.blur(0.1)
                    try:
                        l = 3 if p == 0.5 else 1
                        vppc.show(rd.axvp, color="b", linewidth=l)
                        vspc.show(rd.axvs, color="b", linewidth=l)
                        rhpc.show(rd.axrh, color="b", linewidth=l)
                        prpc.show(rd.axpr, color="b", linewidth=l)

                        #dmout = depthmodel(vppc, vspc, rhpc) #use extract for saveing
                        #dmout.write96('_HerrMet.p%.2f.mod96' % p, overwrite=True)#use extract for saveing
                    except KeyboardInterrupt: raise
                    except Exception as e:
                        print "Error", str(e)
        else:
            print "call option --run to start inversion"

        # ------
        if os.path.exists('_HerrMet.param'):
            p, _ = load_paramfile('_HerrMet.param')
            showvp, showvs, showrh, showpr = True, True, True, True
            if   isinstance(p, Parameterizer_mZVSVPRH): showpr = False
            elif isinstance(p, Parameterizer_mZVSPRRH): showvp = False
            elif isinstance(p, Parameterizer_mZVSPRzRHvp): showvp = showpr = showrh = False
            elif isinstance(p, Parameterizer_mZVSPRzRHz):  showvp = showpr = showrh = False
            else: raise Exception('')

            # inexact
            # rd.plotmodel(alpha=1.0, color="r", marker = "o--", linewidth=3,
            #              showvp=showvp, showvs=showvs, showrh=showrh, showpr=showpr, *p.inv(p.MINF))
            # rd.plotmodel(alpha=1.0, color="r", marker = "o--", linewidth=3,
            #              showvp=showvp, showvs=showvs, showrh=showrh, showpr=showpr, *p.inv(p.MSUP))
            #
            vplow, vphgh, vslow, vshgh, rhlow, rhhgh, prlow, prhgh = p.boundaries()
            vplow.show(rd.axvp, alpha = 1.0, color = "r", marker = "o--", linewidth = 3)
            vphgh.show(rd.axvp, alpha = 1.0, color = "r", marker = "o--", linewidth = 3)
            vslow.show(rd.axvs, alpha = 1.0, color = "r", marker = "o--", linewidth = 3)
            vshgh.show(rd.axvs, alpha = 1.0, color = "r", marker = "o--", linewidth = 3)
            rhlow.show(rd.axrh, alpha = 1.0, color = "r", marker = "o--", linewidth = 3)
            rhhgh.show(rd.axrh, alpha = 1.0, color = "r", marker = "o--", linewidth = 3)
            prlow.show(rd.axpr, alpha = 1.0, color = "r", marker = "o--", linewidth = 3)
            prhgh.show(rd.axpr, alpha = 1.0, color = "r", marker = "o--", linewidth = 3)
            zmax = 1.1 * p.inv(p.MINF)[0][-1]

            if isinstance(p, Parameterizer_mZVSPRzRHvp):
                rd.axpr.plot(
                    p.PRz(np.linspace(0., zmax, 100)),
                    np.linspace(0., zmax, 100), "r--", linewidth = 3)
                textonly(rd.axpr, p.PRzName, loc=4)
                textonly(rd.axrh, p.RHvpName, loc=4)
            elif isinstance(p, Parameterizer_mZVSPRzRHz):
                rd.axpr.plot(
                    p.PRz(np.linspace(0., zmax, 100)),
                    np.linspace(0., zmax, 100), "r--", linewidth = 3)
                rd.axrh.plot(
                    p.RHz(np.linspace(0., zmax, 100)),
                    np.linspace(0., zmax, 100), "r--", linewidth = 3)
                textonly(rd.axpr, p.PRzName, loc=4)
                textonly(rd.axrh, p.RHzName, loc=4)

            rd.set_zlim(np.array([0, zmax]))
        else:
            print "call option --param to see prior depth boundaries"

        # --------------------
        if "m96" in argv.keys():  # plot personal data on top
            for m96 in argv['m96']:
                try:
                    dm = depthmodel_from_mod96(m96)
                    dm.vp.show(rd.axvp, "g", linewidth=3)
                    dm.vs.show(rd.axvs, "g", linewidth=3)
                    dm.rh.show(rd.axrh, "g", linewidth=3)
                    dm.pr().show(rd.axpr, "g", linewidth=3)
                except KeyboardInterrupt: raise
                except Exception as e:
                    print 'could not read or display %s (reason %s)' % (m96, str(e))
        # --------------------
        if "sltz" in argv.keys():  # plot personal data on top (private option)
            dm = depthmodel_from_mod96('/home/max/progdat/CPiS/EarthModel/Soultz.rho.mod')
            dm.vp.show(rd.axvp, "c", linewidth=3)
            dm.vs.show(rd.axvs, "c", linewidth=3)
            dm.rh.show(rd.axrh, "c", linewidth=3)
            dm.pr().show(rd.axpr, "c", linewidth=3)
        # --------------------
        if "ritt" in argv.keys(): # plot personal data on top (private option)
            A = Allfile("/home/max/data/puits/GRT1/GRT1.logsonic")
            A.read()
            rd.axvp.plot(A.data['VP'], A.data['TVD'] / 1000., color="m", alpha=0.4)
            rd.axvs.plot(A.data['VS'], A.data['TVD'] / 1000., color="m", alpha=0.4)
            rd.axpr.plot(A.data['VP'] / A.data['VS'], A.data['TVD'] / 1000., color="m", alpha=0.4)
            dm = depthmodel_from_mod96('/home/max/progdat/CPiS/EarthModel/GRT1.Maurer2016.rho.mod')
            dm.vp.show(rd.axvp, "m", linewidth=3)
            dm.vs.show(rd.axvs, "m", linewidth=3)
            dm.rh.show(rd.axrh, "m", linewidth=3)
            dm.pr().show(rd.axpr, "m", linewidth=3)

        # ------
        if os.path.exists("_HerrMet.target"):
            # plot data on top
            rd.plotdisp(d.waves, d.types, d.modes, d.freqs, d.inv(dobs), dvalues=d.dvalues, alpha=0.8, color="g",linewidth=3)

            if "overdisp" in argv.keys():
                rd.set_vlim((0.5 * d.values.min(), 1.5 * d.values.max()))
                rd.set_plim((0.8 / overfreqs.max(), 1.2 / overfreqs.min()))
            else:
                rd.set_vlim((0.8 * d.values.min(), 1.2 * d.values.max()))
                rd.set_plim((0.8 / d.freqs.max(), 1.2 / d.freqs.min()))
        rd.tick()
        rd.grid()
        chftsz(rd.fig, 16)
        if "png" in argv.keys():
            rd.fig.savefig('_HerrMet.png')
        else:
            showme()

    # # -------------------------------------
    # if "top" in argv.keys() or "pdf" in argv.keys():
    #     # assume called after param, target and run
    #
    #     # ------
    #     if os.path.exists("_HerrMet.target"):
    #         rd = resultdisplay1(targetfile="_HerrMet.target")
    #         d = makedatacoder("./_HerrMet.target", which=Datacoder_log)  # datacoder based on observations
    #         dobs, _ = d.target()
    #     else:
    #         rd = resultdisplay1()
    #         print "call option --target to see the target dispersion curves"
    #
    #     # ------
    #     if os.path.exists('_HerrMet.param'):
    #         p = Parameterizer_mZVSPRRH('_HerrMet.param')
    #         rd.plotmodel(alpha=1.0, color="r", linewidth=3, *p.inv(p.MINF))
    #         rd.plotmodel(alpha=1.0, color="r", linewidth=3, *p.inv(p.MSUP))
    #         zmax = 1.1 * p.inv(p.MINF)[0][-1]
    #         rd.set_zlim(np.array([0, zmax]))
    #     else:
    #         print "call option --param to see prior depth boundaries"
    #     # ------
    #     if os.path.exists('_HerrMet.run'):
    #         if "top" in argv.keys():
    #             top = int(argv['top'][0]) if argv['top'] is not None else default_top
    #             topstep = int(argv['top'][1]) if argv['top'] is not None and len(argv['top']) >= 2 else default_topstep
    #
    #             chainids, weights, llks, ms, ds = readHerrmetout_1('_HerrMet.run', top=top, topstep=topstep)
    #             vmin, vmax, cmap = llks.min(), llks.max(), tej()  # plt.cm.hot#
    #             colors = values2colors(llks, vmin=vmin, vmax=vmax, cmap=cmap)
    #             for i in range(len(llks))[::-1]:
    #                 rd.plotmodel(color=colors[i], alpha=1.0, *ms[i])
    #                 rd.plotdisp(color=colors[i], alpha=1.0, *ds[i])
    #
    #             cb = makecolorbar(vmin=vmin, vmax=vmax, cmap=cmap)
    #             cax = rd.fig.add_axes((0.78, 0.1, 0.005, 0.3))
    #             rd.fig.colorbar(cb, cax=cax, label="log likelyhood")
    #
    #         if "pdf" in argv.keys():
    #             dms, weights = [], []
    #             for chainid, weight, llk, (ztop, vp, vs, rh), (waves, types, modes, freqs, values) in readHerrmetout(
    #                     "_HerrMet.run"):
    #                 dm = depthmodel(depthmodel1D(ztop, vp),
    #                                 depthmodel1D(ztop, vs),
    #                                 depthmodel1D(ztop, rh))
    #                 dms.append(dm)
    #                 weights.append(weight)
    #
    #             for p, (vppc, vspc, rhpc, prpc) in dmstats1(dms, percentiles=[0.16, 0.5, 0.84], Ndepth=100, Nvalue=100,
    #                                                         weights=weights):
    #                 try:
    #                     vppc.show(rd.axvp, color="k", linewidth=2)
    #                     vspc.show(rd.axvs, color="k", linewidth=2)
    #                     rhpc.show(rd.axrh, color="k", linewidth=2)
    #                     prpc.show(rd.axpr, color="k", linewidth=2)
    #
    #                     dmout = depthmodel(vppc, vspc, rhpc)
    #                     dmout.write96('_HerrMet.p%.2f.mod96' % p, overwrite=True)
    #                 except KeyboardInterrupt:
    #                     raise
    #                 except Exception as e:
    #                     print e
    #
    #         rd.tick()
    #
    #     else:
    #         print "call option --run to start inversion"
    #
    #     # --------------------
    #     if "sltz" in argv.keys():  # plot personal data on top
    #         dm = depthmodel_from_mod96('/home/max/progdat/CPiS/EarthModel/Soultz.rho.mod')
    #         dm.vp.show(rd.axvp, "g", linewidth=3)
    #         dm.vs.show(rd.axvs, "g", linewidth=3)
    #         dm.rh.show(rd.axrh, "g", linewidth=3)
    #         dm.pr().show(rd.axpr, "g", linewidth=3)
    #     # -------------------- #plot personal data on top
    #     if "ritt" in argv.keys():
    #         A = Allfile("/home/max/data/puits/GRT1/GRT1.logsonic")
    #         A.read()
    #         rd.axvp.plot(A.data['VP'], A.data['TVD'] / 1000., color="k", alpha=0.5)
    #         rd.axvs.plot(A.data['VS'], A.data['TVD'] / 1000., color="k", alpha=0.5)
    #         rd.axpr.plot(A.data['VP'] / A.data['VS'], A.data['TVD'] / 1000., color="k", alpha=0.5)
    #
    #     # ------
    #     if os.path.exists("_HerrMet.target"):
    #         # plot data on top
    #         rd.plotdisp(d.waves, d.types, d.modes, d.freqs, d.inv(dobs), dvalues=d.dvalues, alpha=0.8, color="k",
    #                     linewidth=2)
    #         rd.set_plim((0.8 / d.freqs.max(), 1.2 / d.freqs.min()))
    #         rd.set_vlim((0.8 * d.values.min(), 1.2 * d.values.max()))
    #         rd.tick()
    #
    #     rd.grid()
    #     if "png" in argv.keys():
    #         rd.fig.savefig('_HerrMet.png')
    #     else:
    #         showme()
