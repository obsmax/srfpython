from srfpython.inversion.metropolis2 import LogUniND
import numpy as np


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
