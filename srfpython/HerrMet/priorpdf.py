from srfpython.inversion.metropolis2 import LogUniND, LogGaussNDCov
import numpy as np


# -------------------------------------
# Prior probability density functions
# -------------------------------------
class DefaultLogRhoM(LogUniND):
    """build the prior pdf, include constraints on the model parameters and eventually on some relations between them
    """
    k, nanbehavior = 1000., 2

    def __init__(self, parameterizer):
        """no more constraints than parameters"""
        LogUniND.__init__(self,
                      vinfs=parameterizer.MINF,
                      vsups=parameterizer.MSUP,
                      k=self.k,
                      nanbehavior=self.nanbehavior)
        self.p = parameterizer

    @staticmethod
    def header(dvp=None, dvs=None, drh=None, dpr=None):
        """ a method to write the apropriate header lines in the parameter file, see subclasses"""
        assert dvp is None
        assert dvs is None
        assert drh is None
        assert dpr is None
        return ""


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

    def __call__(self, m):
        _, _, VS, _ = self.p.inv(m)
        DVS = VS[1:] - VS[:-1]
        extended_m = np.concatenate((m, DVS))
        return LogUniND.__call__(self, extended_m)

    @staticmethod
    def header(dvp=None, dvs=None, drh=None, dpr=None):
        assert dvp is None
        assert drh is None
        assert dpr is None
        header = '#met PRIORTYPE = "DVS"\n'
        header += '#met DVSMIN = {}\n'.format(dvs[0])
        header += '#met DVSMAX = {}\n'.format(dvs[1])
        return header


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

    def __call__(self, m):
        _, VP, VS, RH = self.p.inv(m)
        DVP = VS[1:] - VS[:-1]
        DVS = VP[1:] - VP[:-1]
        DRH = RH[1:] - RH[:-1]
        extended_m = np.concatenate((m, DVP, DVS, DRH))
        return LogUniND.__call__(self, extended_m)

    @staticmethod
    def header(dvp=None, dvs=None, drh=None, dpr=None):
        assert dpr is None
        header  = '#met PRIORTYPE = "DVPDVSDRH"\n'
        header += '#met DVPMIN = {}\n'.format(dvp[0])
        header += '#met DVPMAX = {}\n'.format(dvp[1])
        header += '#met DVSMIN = {}\n'.format(dvs[0])
        header += '#met DVSMAX = {}\n'.format(dvs[1])
        header += '#met DRHMIN = {}\n'.format(drh[0])
        header += '#met DRHMAX = {}\n'.format(drh[1])
        return header


class LogRhoM_DVPDVSDRHDPR(DefaultLogRhoM):
    """add new constraitns on the vs offsets on interfaces"""

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

    def __call__(self, m):
        _, VP, VS, RH = self.p.inv(m)
        DVP = VS[1:] - VS[:-1]
        DVS = VP[1:] - VP[:-1]
        DRH = RH[1:] - RH[:-1]
        PR = VP / VS
        DPR = PR[1:] - PR[:-1]
        extended_m = np.concatenate((m, DVP, DVS, DRH, DPR))
        return LogUniND.__call__(self, extended_m)

    @staticmethod
    def header(dvp=None, dvs=None, drh=None, dpr=None):

        header  = '#met PRIORTYPE = "DVPDVSDRHDPR"\n'
        header += '#met DVPMIN = {}\n'.format(dvp[0])
        header += '#met DVPMAX = {}\n'.format(dvp[1])
        header += '#met DVSMIN = {}\n'.format(dvs[0])
        header += '#met DVSMAX = {}\n'.format(dvs[1])
        header += '#met DRHMIN = {}\n'.format(drh[0])
        header += '#met DRHMAX = {}\n'.format(drh[1])
        header += '#met DPRMIN = {}\n'.format(dpr[0])
        header += '#met DPRMAX = {}\n'.format(dpr[1])
        return header


class LogRhoM_TIKHONOV(DefaultLogRhoM):
    """add new constraitns on the vs offsets on interfaces"""

    def __init__(self, parameterizer, alpha_tikhonov: float):
        DefaultLogRhoM.__init__(self, parameterizer=parameterizer)
        self.alpha_tikhonov = alpha_tikhonov

        # compute the range of each parameters to normalize the gradients
        (vplow, vphgh,
         vslow, vshgh,
         rhlow, rhhgh,
         prlow, prhgh) = self.p.boundaries()

        self.vs_range = vshgh.values.max() - vslow.values.min()
        self.vp_range = vphgh.values.max() - vplow.values.min()
        self.rh_range = rhhgh.values.max() - rhlow.values.min()
        self.pr_range = prhgh.values.max() - prlow.values.min()

    def __call__(self, m):

        log_rhom_1 = DefaultLogRhoM.__call__(self, model=m)

        ztop, vp, vs, rh = self.p.inv(m)
        pr = vp / vs

        # top of half space
        depthmax = ztop[-1]

        # layer thicknesses
        dz = ztop[1:] - ztop[:-1]

        # vertical gradients in each layer except half space
        dvp_o_dz = (vs[1:] - vs[:-1]) / dz
        dvs_o_dz = (vp[1:] - vp[:-1]) / dz
        drh_o_dz = (rh[1:] - rh[:-1]) / dz
        dpr_o_dz = (pr[1:] - pr[:-1]) / dz

        # scale the gradients
        dvs_o_dz *= depthmax / self.vs_range
        dvp_o_dz *= depthmax / self.vp_range
        drh_o_dz *= depthmax / self.rh_range
        dpr_o_dz *= depthmax / self.pr_range        #

        # gradients norms
        log_rhom_2 = -0.5 * (dvs_o_dz ** 2.0).sum()
        log_rhom_3 = -0.5 * (dvp_o_dz ** 2.0).sum()
        log_rhom_4 = -0.5 * (drh_o_dz ** 2.0).sum()
        log_rhom_5 = -0.5 * (dpr_o_dz ** 2.0).sum()

        return log_rhom_1 +\
            self.alpha_tikhonov * (log_rhom_2 + log_rhom_3 + log_rhom_4 + log_rhom_5)

    @staticmethod
    def header(dvp=None, dvs=None, drh=None, dpr=None):

        header  = '#met PRIORTYPE = "DVPDVSDRHDPR"\n'
        header += '#met DVPMIN = {}\n'.format(dvp[0])
        header += '#met DVPMAX = {}\n'.format(dvp[1])
        header += '#met DVSMIN = {}\n'.format(dvs[0])
        header += '#met DVSMAX = {}\n'.format(dvs[1])
        header += '#met DRHMIN = {}\n'.format(drh[0])
        header += '#met DRHMAX = {}\n'.format(drh[1])
        header += '#met DPRMIN = {}\n'.format(dpr[0])
        header += '#met DPRMAX = {}\n'.format(dpr[1])
        return header