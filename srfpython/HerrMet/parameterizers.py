import numpy as np
from srfpython.depthdisp.depthmodels import depthmodel1D
from srfpython.utils import string2func
import os
r43 = np.sqrt(4./3.)

"""
Parameterizer : object to convert an array of parameters into smth that can be understood by Herrmann.dispersion
Theory        : object to convert a model array into a data array
Datacoder     : object to convert output from Herrmann.dispersion to an array of data
                          a model array (m)                        
  |                          ^    |
  |                          |    |
  |      mod96file ----->  parameterizer   --------> m_apr, CM
  |      (apriori)           |    |
  |                          |    v
  |                        depth model 
  |                     (ztop, vp, vs, rh)
  |                            |
 theory                 Herrmann.dispersion
 (forward problem)             |
  |                            v
  |                       dispersion data 
  |                      (waves, types, modes, freqs, values, (dvalues))
  |                          ^    |
  |                          |    |
  |     surf96file ----->   datacoder      --------> d_obs, CD => logRHOD
  |      (target)            |    |
  v                          |    v
                          a data array (d)
"""


def check_parameter_file(A):
    """
    :param A: a read AsciiFile object corresponding to a parameter file
    :return:
    """
    metakeys = A.metadata.keys()
    if not "NLAYER" in metakeys:
        raise ValueError('NLAYER not found in file metadata')
    if not "TYPE" in metakeys:
        raise ValueError('TYPE not found in file metadata')
    if not len(A.data['KEY']) == len(np.unique(A.data['KEY'])):
        raise ValueError('there are repeated entries in column KEYS')
    if np.any(A.data['VINF'] > A.data['VSUP']):
        raise ValueError('VSUP cannot be lower than VINF')

    # TODO add more verifications in common to all parameterizers here


# -------------------------------------
class Parameterizer(object):
    """the parameterizer object links the model array (m) to a set of variables that can
       be understood by the theory function (ZTOP, VP, VS, RH)
       this is the default class, methods must be overwritten by subclasses"""

    # ------------------
    def __init__(self, A):
        """
        :param A: AsciiFile with parameterization

        please set up attributes :
        self.NLAYER : int, the number of layers including half space
        self.MDEFAULT : array, all the parameters including those that should not be inverted (not-masked)
        self.I : boolean array, mask used to determine the parameters to be inverted, the others will
                 keep constant and equal to self.MDEFAULT
        self.MMEAN = array, the mean model, masked by self.I
        self.MINF = array, lower boundary for each parameter, masked by self.I
        self.MSUP = array, upper boundary for each parameter, masked by self.I
        self.MSTD = array, markov proposal for each parameter, masked by self.I
        """
        check_parameter_file(A)
        raise Exception('Never used : please custom subclasses')

    # ------------------
    def boundaries(self):
        """
        a method to determine the lowest and highest models, no arguments
        based on columns VINF and VSUP in the parameterization file
        warning : the lower model is not necessarily obtained when all parameters reach their lowest boundary (VINF)
        please return lowest and highest depthmodels stored as depthmodel1D
        vplow, vphgh,
        vslow, vshgh,
        rhlow, rhhgh,
        prlow, prhgh
        """

        # default behavior, to be customized if needed
        Ztopsup, VPlow, VSlow, RHlow = self.inv(self.MINF)
        Ztopinf, VPhgh, VShgh, RHhgh = self.inv(self.MSUP)
        z = np.sort(np.unique(np.concatenate((Ztopinf, Ztopsup))))

        # --------------------
        def f(Zinf, Zsup, V, which):
            v1 = depthmodel1D(Zinf, V).interp(z, interpmethod="stairs")
            v2 = depthmodel1D(Zsup, V).interp(z, interpmethod="stairs")
            return depthmodel1D(z, which(np.concatenate(([v1], [v2]), axis=0), axis=0))

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
        """a method to provide the name of the inverted parameters

        please return the masked version of the array, to return only the keys that
        are actually inverted (use the boolean array self.I)

        e.g. keys = np.array(['Z1', 'VP0', 'VP1', 'VS0', 'VS1', 'RH0', 'RH1'])[self.I]

        :return: array or list of parameters
        """
        raise Exception('never used, please custom subclasses')

    # ------------------
    def __call__(self, ZTOP, VP, VS, RH):
        """the method that is called to convert
           a real depth model (ZTOP, VP, VS, RH)
           into a model array (m = array of parameters).
           Called inside Theory.__call__ just before dispersion
           please keep consistent with self.inv
           use self.keys to get the corresponding parameter names

           :return: array of parameters m (corresponding to self.I)
        """
        raise Exception('never used, please custom subclasses')

    # ------------------
    def inv(self, m):
        """the method that is called to convert model array
           into a true depth model that can be understood by srfdis17.
           Called for plotting (i.e. convert a model back to something more physical)
           please keep consistent with self.__call__
           :return: 4 arrays ZTOP, VP, VS, RH
        """
        raise Exception('never used, please custom subclasses')

    # ------------------
    def prior(self, ZTOP_apr, VP_apr, VS_apr, RH_apr):
        """
        to be used for linearized inversion (HerrLin)
        :return: prior model (mapr) and covariance matrix (CMinv)
                 in agreement with the parameterizer type
        """
        raise Exception('never used, please custom subclasses')


# -------------------------------------
class Parameterizer_mZVSPRRH(Parameterizer):

    def __init__(self, A):
        """see Parameterizer"""
        check_parameter_file(A)
        assert A.metadata['TYPE'] == "mZVSPRRH"
        assert np.all(A.data['VINF'] <= A.data['VSUP'])
        if np.all(A.data['VINF'] == A.data['VSUP']):
            print "Warning : all parameters are locked"
        # --------
        self.NLAYER = A.metadata['NLAYER']
        self.I = A.data['VINF'] < A.data['VSUP'] #index of dimension used, other parameters are kept constant
        # --------
        self.MDEFAULT = 0.5 * (A.data['VINF'] + A.data['VSUP']) #used for dimensions that are not part of the model
        self.MMEAN = self.MDEFAULT[self.I]
        self.MINF  = A.data['VINF'][self.I]
        self.MSUP  = A.data['VSUP'][self.I]
        self.MSTD  = 0.5 * (A.data['VSUP'] - A.data['VINF'])[self.I]
        # --------

    # ------------------
    def boundaries(self):
        """see Parameterizer"""
        Ztopsup, VPlow, VSlow, RHlow = self.inv(self.MINF)
        Ztopinf, VPhgh, VShgh, RHhgh = self.inv(self.MSUP)
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
        """see Parameterizer"""
        k = ["-Z%d" % i for i in xrange(1, self.NLAYER)] + \
             ['VS%d' % i for i in xrange(self.NLAYER)] + \
             ['PR%d' % i for i in xrange(self.NLAYER)] + \
             ['RH%d' % i for i in xrange(self.NLAYER)]
        return np.asarray(k)[self.I]

    # ------------------
    def __call__(self, ZTOP, VP, VS, RH):
        """see Parameterizer"""
        m = np.concatenate((-1.0 * ZTOP[1:], VS, VP / VS, RH))[self.I]
        return m

    # ------------------
    def inv(self, m):
        """see Parameterizer"""
        M = self.MDEFAULT.copy()
        M[self.I] = m #overwrites default values with the one provided in m

        nlayer = self.NLAYER#(len(M) + 1) / 4
        ZTOP = np.concatenate(([0.], -1.0 * M[: 1 * nlayer - 1]))
        VS = M[1 * nlayer - 1: 2 * nlayer - 1]
        PR = M[2 * nlayer - 1: 3 * nlayer - 1]
        RH = M[3 * nlayer - 1: 4 * nlayer - 1]
        VP = PR * VS

        return ZTOP, VP, VS, RH

    # ------------------
    def prior(self, ZTOP_apr, VP_apr, VS_apr, RH_apr,
              **kwargs):
        """see Parameterizer"""
        mapr = self.__call__(ZTOP_apr, VP_apr, VS_apr, RH_apr)
        raise NotImplementedError('')


# -------------------------------------
class Parameterizer_mZVSVPRH(Parameterizer):
    """see Parameterizer_mZVSPRRH for doc"""
    def __init__(self, A):
        check_parameter_file(A)
        assert A.metadata['TYPE'] == "mZVSVPRH"
        assert np.all(A.data['VINF'] <= A.data['VSUP'])
        if np.all(A.data['VINF'] == A.data['VSUP']):
            print "Warning : all parameters are locked"
        # --------
        self.NLAYER = A.metadata['NLAYER']
        self.I = A.data['VINF'] < A.data['VSUP']
        # --------
        self.MDEFAULT = 0.5 * (A.data['VINF'] + A.data['VSUP'])      # value to use if not a parameter
        self.MMEAN = self.MDEFAULT[self.I]                     #
        self.MINF = A.data['VINF'][self.I]                        # lower boundary for each parameter
        self.MSUP = A.data['VSUP'][self.I]                        # upper boundary for each parameter
        self.MSTD = 0.5 * (A.data['VSUP'] - A.data['VINF'])[self.I]  # markov proposal for each parameter
        # --------

    # ------------------
    def boundaries(self):
        """see Parameterizer"""
        Ztopsup, VPlow, VSlow, RHlow = self.inv(self.MINF)
        Ztopinf, VPhgh, VShgh, RHhgh = self.inv(self.MSUP)
        z = np.sort(np.unique(np.concatenate((Ztopinf, Ztopsup))))

        # --------------------
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
        """see Parameterizer"""
        k = ["-Z%d" % i for i in xrange(1, self.NLAYER)] + \
            ['VS%d' % i for i in xrange(self.NLAYER)] + \
            ['VP%d' % i for i in xrange(self.NLAYER)] + \
            ['RH%d' % i for i in xrange(self.NLAYER)]
        return np.asarray(k)[self.I]

    # ------------------
    def __call__(self, ZTOP, VP, VS, RH):
        """see Parameterizer"""
        m = np.concatenate((-1.0 * ZTOP[1:], VS, VP, RH))[self.I]
        return m

    # ------------------
    def inv(self, m):
        """see Parameterizer"""
        M = self.MDEFAULT.copy()
        M[self.I] = m  # overwrites default values with the one provided in m

        nlayer = self.NLAYER  # (len(M) + 1) / 4
        ZTOP = np.concatenate(([0.], -1.0 * M[: 1 * nlayer - 1]))
        VS = M[1 * nlayer - 1: 2 * nlayer - 1]
        VP = M[2 * nlayer - 1: 3 * nlayer - 1]
        RH = M[3 * nlayer - 1: 4 * nlayer - 1]

        return ZTOP, VP, VS, RH

    # ------------------
    def prior(self, ZTOP_apr, VP_apr, VS_apr, RH_apr,
              **kwargs):
        """see Parameterizer"""
        mapr = self.__call__(ZTOP_apr, VP_apr, VS_apr, RH_apr)

        raise NotImplementedError('')


# -------------------------------------
class Parameterizer_mZVSPRzRHvp(Parameterizer):
    """see Parameterizer_mZVSPRRH for doc"""

    def __init__(self, A):
        """see Parameterizer"""
        check_parameter_file(A)
        assert A.metadata['TYPE'] == "mZVSPRzRHvp"
        assert np.all(A.data['VINF'] <= A.data['VSUP'])
        if np.all(A.data['VINF'] == A.data['VSUP']):
            print "Warning : all parameters are locked"
        # --------
        self.NLAYER = A.metadata['NLAYER']
        self.I = A.data['VINF'] < A.data['VSUP']
        # --------
        self.MDEFAULT = 0.5 * (A.data['VINF'] + A.data['VSUP'])
        self.MMEAN = self.MDEFAULT[self.I]
        self.MINF = A.data['VINF'][self.I]
        self.MSUP = A.data['VSUP'][self.I]
        self.MSTD = 0.5 * (A.data['VSUP'] - A.data['VINF'])[self.I]
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

    # ------------------
    def boundaries(self):
        """see Parameterizer"""
        print "boundaries method not implemented for %s, using default one" % self.__class__.__name__
        return Parameterizer.boundaries(self)

    # ------------------
    def keys(self):
        """see Parameterizer"""
        k = ["-Z%d" % i for i in xrange(1, self.NLAYER)] + \
            ['VS%d' % i for i in xrange(self.NLAYER)]
        return np.asarray(k)[self.I]

    # ------------------
    def __call__(self, ZTOP, VP, VS, RH):
        """see Parameterizer"""
        m = np.concatenate((-1.0 * ZTOP[1:], VS))[self.I]
        return m

    # ------------------
    def inv(self, m):
        """see Parameterizer"""
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

    # ------------------
    def prior(self, ZTOP_apr, VP_apr, VS_apr, RH_apr, **kwargs):
        """see Parameterizer"""
        raise NotImplementedError('')


# -------------------------------------
class Parameterizer_mZVSPRzRHz(Parameterizer):
    """see Parameterizer_mZVSPRRH for doc"""

    def __init__(self, A):
        """see Parameterizer"""
        check_parameter_file(A)
        assert A.metadata['TYPE'] == "mZVSPRzRHz"
        assert np.all(A.data['VINF'] <= A.data['VSUP'])
        if np.all(A.data['VINF'] == A.data['VSUP']):
            print "Warning : all parameters are locked"
        # --------
        self.NLAYER = A.metadata['NLAYER']
        self.I = A.data['VINF'] < A.data['VSUP']
        # --------
        self.MDEFAULT = 0.5 * (A.data['VINF'] + A.data['VSUP'])
        self.MMEAN = self.MDEFAULT[self.I]
        self.MINF = A.data['VINF'][self.I]
        self.MSUP = A.data['VSUP'][self.I]
        self.MSTD = 0.5 * (A.data['VSUP'] - A.data['VINF'])[self.I]
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

    # ------------------
    def boundaries(self):
        print "boundaries method not implemented for %s, using default one" % self.__class__.__name__
        return Parameterizer.boundaries(self)

    # ------------------
    def keys(self):
        """see Parameterizer"""
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
        """see Parameterizer"""
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

    # ------------------
    def prior(self, ZTOP_apr, VP_apr, VS_apr, RH_apr, **kwargs):
        """see Parameterizer"""
        raise NotImplementedError('')

