from __future__ import print_function

#from labex.inversion.linearized6 import *
#from labex.dispersion.depthmodels import *
#from labex.dispersion.displaws import *
#from labex.dispersion.Herrmann0.dispersion import dispersion as srfdis17, dispersion_1 as srfdis17_1
#from labex.graphictools.gutils import plt, grapher, gcf, gca, showme, pause, logtick, Ntick, value2color
#from labex.graphictools.cmaps import  tej
import numpy as np
r43 = np.sqrt(4. / 3.)



"""
                          a model array (m)                        
  |                          ^    |
  |                          |    |
  |      mod96file ----->  parameterizer   --------> m_apr, CM
  |      (apriori)           |    |
  |                          |    v
  |                        depth model 
  |                     (ztop, vp, vs, rh)
  |                            |
 theory                     srfdis17 (modified Herrmann)
 (forward problem)             |
  |                            v
  |                       dispersion data 
  |                      (waves, types, modes, freqs, values, (dvalues))
  |                          ^    |
  |                          |    |
  |     surf96file ----->   datacoder      --------> d_obs, CD
  |      (target)            |    |
  v                          |    v
                          a data array (d)
"""

# ------------------
class Parameterizer(object):
    def __init__(self, ZTOP, VP, VS, RH):
        "init with the starting model and uncertainty"
        self.ZTOP = ZTOP
        self.VP   = VP
        self.VS   = VS
        self.RH   = RH

    # ------------------
    def getrange(self, hmin = 0.001, hmax = 1.0, prmin = r43, prmax = 2.5, vsmin = 0.08, vsmax = 4.0, rhmin = 2.0, rhmax = 3.0):
        "determines the lower and upper bounds in each dimension"
        assert  0.001 <= hmin  < hmax
        assert  r43   <= prmin < prmax
        assert  0.08  <= vsmin < vsmax
        assert  1.0   <= rhmin < rhmax
        nlayer = len(self.VS)
        o = np.ones(nlayer, float)
        ztopinf = np.concatenate(([0.], (hmin * o)[:-1].cumsum()))
        ztopsup = np.concatenate(([0.], (hmax * o)[:-1].cumsum()))
        minf = self(ztopinf, prmin * vsmin * o, vsmin * o, rhmin * o)
        msup = self(ztopsup, prmax * vsmax * o, vsmax * o, rhmax * o)
        return minf, msup        
#------------------
class Parameterizer_logH_logMU_logK_logRH(Parameterizer):
    def __call__(self, ZTOP, VP, VS, RH):
        H  = ZTOP[1:] - ZTOP[:-1]
        MU = RH * VS ** 2.
        K  = RH * VP ** 2. - 4. * MU / 3.
        return np.log(np.concatenate((H, MU, K, RH)))
    #------------------
    def inv(self, m):
        assert not (len(m) + 1) % 4
        nlayer = (len(m) + 1) / 4
        logH   = m[              : 1 * nlayer -1]
        logMU  = m[1 * nlayer - 1: 2 * nlayer -1]
        logK   = m[2 * nlayer - 1: 3 * nlayer -1]
        logRH  = m[3 * nlayer - 1:]

        H    = np.exp(logH)
        MU   = np.exp(logMU)
        K    = np.exp(logK)
        RH   = np.exp(logRH)
        ZTOP = np.concatenate(([0.], H.cumsum()))
        VP   = np.sqrt((K + 4. * MU / 3.) / RH)
        VS   = np.sqrt(MU / RH)
        return ZTOP, VP, VS, RH
    #------------------
    def startmodel(self, sigH = 0.01, sigVP = 0.1, sigVS = 0.1, sigRH = 0.1, smoothMU = 0., smoothK = 0., smoothRH = 0., lockbot = False):#sigz = 0.01, sigvs = 0.1, sigpr = 0.0001, sigrh = 0.1, smoothmu = 0., smoothK = 0., smoothrh = 0.):
        mapr  = self.__call__(self.ZTOP, self.VP, self.VS, self.RH)
        H = self.ZTOP[1:] - self.ZTOP[:-1]
        siglogH = sigH / abs(H)
        siglogRH = sigRH / self.RH
        siglogMU = siglogRH + 2. * sigVS / self.VS
        toto = abs(self.VP ** 2. - (4. / 3.) * self.VS ** 2.)
        siglogK  = siglogRH + 2. * sigVP * self.VP / toto  + (8. / 3.) * self.VS * sigVS / toto

        if lockbot:
            assert len(sigVS) > 1
            #cannot lock sigH
            siglogMU[-1] = siglogK[-1] = siglogRH[-1] = 1.e-20

        CMinv = np.concatenate((siglogH, siglogMU, siglogK, siglogRH)) ** -2.

        if np.any([smoothMU, smoothK, smoothRH]):
            N = len(self.VS)
            rho = np.eye(len(mapr))
            r    = np.arange(N)
            I, J = np.meshgrid(r, r)
            D = abs(I - J)
            if smoothMU:
                rho[N-1:2*N-1, N-1:2*N-1] = np.exp(-D / smoothMU) #smoothing distance in cell number
            if smoothK:
                rho[2*N-1:3*N-1, 2*N-1:3*N-1] = np.exp(-D / smoothK) #smoothing distance in cell number
            if smoothRH:
                rho[3*N-1:4*N-1, 3*N-1:4*N-1] = np.exp(-D / smoothRH) #smoothing distance in cell number

            if False:
                plt.figure()
                gca().matrixview(rho)
                showme()


            VX, VY = np.meshgrid(1. / CMinv, 1. / CMinv)
            CM = rho * np.sqrt(VX * VY)
            CMinv = np.linalg.inv(CM)

        return mapr, CMinv
#------------------
class Parameterizer_logH_logVS_logPR_logRH(Parameterizer):
    """ """
    def __call__(self, ZTOP, VP, VS, RH):
        H  = ZTOP[1:] - ZTOP[:-1]

        m = np.concatenate((H - (0.001 - 0.0001), VS - (0.08 - 0.001), VP/VS - (r43 - 1.e-3), RH - (1. - 0.001))) #must be consistent with inv
        I = (m > 0.)
        m[I] = np.log(m[I])
        m[~I] = np.nan
        return m
    #------------------
    def inv(self, m):
        assert not (len(m) + 1) % 4
        nlayer = (len(m) + 1) / 4
        logH   = m[              : 1 * nlayer -1]
        logVS  = m[1 * nlayer - 1: 2 * nlayer -1]
        logPR  = m[2 * nlayer - 1: 3 * nlayer -1]
        logRH  = m[3 * nlayer - 1:]

        H    = np.exp(logH)  + (0.001 - 0.0001) #must be consistent with inv
        VS   = np.exp(logVS) + (0.08 - 0.001)
        PR   = np.exp(logPR) + (r43 - 1.e-3)
        RH   = np.exp(logRH) + (1. - 0.001)

        ZTOP = np.concatenate(([0.], H.cumsum()))
        VP   = VS * PR

        return ZTOP, VP, VS, RH

    #------------------
    def startmodel(self, sigH = 0.01, sigVS = 0.1, sigPR = 0.1, sigRH = 0.1, smoothVS = 0., smoothPR = 0., smoothRH = 0., lockbot = False):#sigz = 0.01, sigvs = 0.1, sigpr = 0.0001, sigrh = 0.1, smoothmu = 0., smoothK = 0., smoothrh = 0.):
        mapr  = self.__call__(self.ZTOP, self.VP, self.VS, self.RH)
        H = self.ZTOP[1:] - self.ZTOP[:-1]
        siglogH = sigH / abs(H)
        siglogVS = sigVS / self.VS
        siglogPR = sigPR / (self.VP / self.VS)
        siglogRH = sigRH / abs(self.RH)

        if lockbot:
            #cannot lock sigH   
            assert len(siglogVS) > 1
            siglogVS[-1] = siglogPR[-1] = siglogRH[-1] = 1.e-20

        CMinv = np.concatenate((siglogH, siglogVS, siglogPR, siglogRH)) ** -2.

        if np.any([smoothVS, smoothPR, smoothRH]):
            N = len(self.VS)
            rho = np.eye(len(mapr))
            r    = np.arange(N)
            I, J = np.meshgrid(r, r)
            D = abs(I - J)
            if smoothVS:
                rho[N-1:2*N-1, N-1:2*N-1] = np.exp(-D / smoothVS) #smoothing distance in cell number
            if smoothPR:
                rho[2*N-1:3*N-1, 2*N-1:3*N-1] = np.exp(-D / smoothPR) #smoothing distance in cell number
            if smoothRH:
                rho[3*N-1:4*N-1, 3*N-1:4*N-1] = np.exp(-D / smoothRH) #smoothing distance in cell number

            if False:
                plt.figure()
                gca().matrixview(rho)
                showme()


            VX, VY = np.meshgrid(1. / CMinv, 1. / CMinv)
            CM = rho * np.sqrt(VX * VY)
            CMinv = np.linalg.inv(CM)

        return mapr, CMinv
#------------------
class Parameterizer_H_VS_PR_RH(Parameterizer):
    """
    """
    def __call__(self, ZTOP, VP, VS, RH):
        H  = ZTOP[1:] - ZTOP[:-1]
        return np.concatenate((H, VS, VP/VS, RH))
    #------------------
    def inv(self, m):
        assert not (len(m) + 1) % 4
        nlayer = (len(m) + 1) / 4
        H   = np.asarray(m[              : 1 * nlayer -1])
        VS  = np.asarray(m[1 * nlayer - 1: 2 * nlayer -1])
        PR  = np.asarray(m[2 * nlayer - 1: 3 * nlayer -1])
        RH  = np.asarray(m[3 * nlayer - 1:])

        ZTOP = np.concatenate(([0.], H.cumsum()))
        VP   = VS * PR

        return ZTOP, VP, VS, RH
    #------------------
    def startmodel(self, sigH = 0.01, sigVS = 0.1, sigPR = 0.1, sigRH = 0.1, smoothVS = 0., smoothPR = 0., smoothRH = 0., lockbot = False):#sigz = 0.01, sigvs = 0.1, sigpr = 0.0001, sigrh = 0.1, smoothmu = 0., smoothK = 0., smoothrh = 0.):
        mapr  = self.__call__(self.ZTOP, self.VP, self.VS, self.RH)
        H = self.ZTOP[1:] - self.ZTOP[:-1]


        sigH  = sigH  * np.ones_like(self.ZTOP[:-1])
        sigVS = sigVS * np.ones_like(self.ZTOP)
        sigPR = sigPR * np.ones_like(self.ZTOP)
        sigRH = sigRH * np.ones_like(self.ZTOP)
        if lockbot:
            #cannot lock sigH   
            assert len(sigVS) > 1
            sigVS[-1] = sigPR[-1] = sigRH[-1] = 1.e-20

        CMinv = np.concatenate((sigH, sigVS, sigPR, sigRH)) ** -2.

        if np.any([smoothVS, smoothPR, smoothRH]):
            N = len(self.VS)
            rho = np.eye(len(mapr))
            r    = np.arange(N)
            I, J = np.meshgrid(r, r)
            D = abs(I - J)
            if smoothVS:
                rho[N-1:2*N-1, N-1:2*N-1] = np.exp(-D / smoothVS) #smoothing distance in cell number
            if smoothPR:
                rho[2*N-1:3*N-1, 2*N-1:3*N-1] = np.exp(-D / smoothPR) #smoothing distance in cell number
            if smoothRH:
                rho[3*N-1:4*N-1, 3*N-1:4*N-1] = np.exp(-D / smoothRH) #smoothing distance in cell number

            if False:
                plt.figure()
                gca().matrixview(rho)
                showme()


            VX, VY = np.meshgrid(1. / CMinv, 1. / CMinv)
            CM = rho * np.sqrt(VX * VY)
            CMinv = np.linalg.inv(CM)

        return mapr, CMinv
#------------------
class Parameterizer_logH_logVS(Parameterizer):
    """
    """
    def __call__(self, ZTOP, VP, VS, RH):
        H  = ZTOP[1:] - ZTOP[:-1]
        return np.log(np.concatenate((H, VS)))
    #------------------
    def inv(self, m):
        assert not (len(m) + 1) % 2
        nlayer = (len(m) + 1) / 2
        logH   = m[              : 1 * nlayer -1]
        logVS  = m[1 * nlayer - 1: 2 * nlayer -1]

        H    = np.exp(logH)
        VS   = np.exp(logVS)

        ZTOP = np.concatenate(([0.], H.cumsum()))
        VP   = VS * (self.VP / self.VS) #keep PR constant
        RH   = self.RH

        return ZTOP, VP, VS, RH
    #------------------
    def startmodel(self, sigH = 0.01, sigVS = 0.1, smoothVS = 0., lockbot = False):
        mapr  = self.__call__(self.ZTOP, self.VP, self.VS, self.RH)
        H = self.ZTOP[1:] - self.ZTOP[:-1]
        siglogH = sigH / abs(H)
        siglogVS = sigVS / self.VS

        if lockbot:
            #cannot lock sigH   
            assert len(siglogVS) > 1
            siglogVS[-1] = 1.e-20

        CMinv = np.concatenate((siglogH, siglogVS)) ** -2.

        if smoothVS:
            N = len(self.VS)
            rho = np.eye(len(mapr))
            r    = np.arange(N)
            I, J = np.meshgrid(r, r)
            D = abs(I - J)
            rho[N-1:2*N-1, N-1:2*N-1] = np.exp(-D / smoothVS) #smoothing distance in cell number
            if False:
                plt.figure()
                gca().matrixview(rho)
                showme()
            VX, VY = np.meshgrid(1. / CMinv, 1. / CMinv)
            CM = rho * np.sqrt(VX * VY)
            CMinv = np.linalg.inv(CM)

        return mapr, CMinv
#------------------
class Parameterizer_logVS(Parameterizer):
    """
    """
    def __call__(self, ZTOP, VP, VS, RH):
        return np.log(VS)
    #------------------
    def inv(self, m):
        assert len(m) == len(self.VS)
        logVS  = m
        VS   = np.exp(logVS)
        VP   = VS * (self.VP / self.VS) #keep PR constant

        return self.ZTOP, VP, VS, self.RH
    #------------------
    def startmodel(self, sigVS = 0.1, smoothVS = 0., lockbot = False):
        mapr  = self.__call__(self.ZTOP, self.VP, self.VS, self.RH)
        siglogVS = sigVS / self.VS

        if lockbot:
            #cannot lock sigH   
            assert len(siglogVS) > 1
            siglogVS[-1] = 1.e-20

        CMinv = siglogVS ** -2.

        if smoothVS:
            N = len(self.VS)
            rho = np.eye(len(mapr))
            r    = np.arange(N)
            I, J = np.meshgrid(r, r)
            D = abs(I - J)
            rho = np.exp(-D / smoothVS) #smoothing distance in cell number
            if False:
                plt.figure()
                gca().matrixview(rho)
                showme()
            VX, VY = np.meshgrid(1. / CMinv, 1. / CMinv)
            CM = rho * np.sqrt(VX * VY)
            CMinv = np.linalg.inv(CM)

        return mapr, CMinv
#------------------
class Parameterizer_logH(Parameterizer):
    """
    """
    def __call__(self, ZTOP, VP, VS, RH):
        H  = ZTOP[1:] - ZTOP[:-1]
        return np.log(H)
    #------------------
    def inv(self, m):
        nlayer = (len(m) + 1)
        logH   = m

        H    = np.exp(logH)
        ZTOP = np.concatenate(([0.], H.cumsum()))
        
        VS   = self.VS
        VP   = self.VP#VS * (self.VP / self.VS) #keep PR constant
        RH   = self.RH

        return ZTOP, VP, VS, RH
    #------------------
    def startmodel(self, sigH = 0.1, lockbot = False):
        mapr  = self.__call__(self.ZTOP, self.VP, self.VS, self.RH)
        H = self.ZTOP[1:] - self.ZTOP[:-1]
        siglogH = sigH / abs(H)

        if lockbot:
            raise Exception('cannot lock ZTOP[-1]')

        CMinv = siglogH ** -2.

        return mapr, CMinv
#------------------
class Parameterizer_logH_logVS_PRlaw_RHlaw(Parameterizer):
    """logH and logVS are parameters,
       VP is recovered as a function of VS using VP/VS = f(depth)
       RH is recovered as a function of VP using RHO = f(VP)
    """
    def __init__(self, ZTOP, VP, VS, RH, prlaw = GRT1pr, rhlaw = gardner74):
        "ignore VP RH"
        
        del VP, RH
        ZMID = np.concatenate((0.5 * (ZTOP[1:] + ZTOP[:-1]), [1.1 * ZTOP[-1]]))

        VP = VS * prlaw(ZMID)
        RH = rhlaw(VP)
        Parameterizer.__init__(self, ZTOP, VP, VS, RH)
        self.prlaw = prlaw
        self.rhlaw = rhlaw
    #------------------
    def __call__(self, ZTOP, VP, VS, RH):
        "ignore VP RH"
        H  = ZTOP[1:] - ZTOP[:-1]
        return np.concatenate((np.log(H), np.log(VS)))
    #------------------
    def inv(self, m):
        assert not (len(m) + 1) % 2
        nlayer = (len(m) + 1) / 2
        logH   = m[              : 1 * nlayer -1]
        logVS  = m[1 * nlayer - 1: 2 * nlayer -1]

        H    = np.exp(logH)
        VS   = np.exp(logVS)
        ZTOP = np.concatenate(([0.], H.cumsum()))
        ZMID = np.concatenate((0.5 * (ZTOP[1:] + ZTOP[:-1]), [1.1 * ZTOP[-1]]))

        VP = VS * self.prlaw(ZMID)
        RH = self.rhlaw(VP)

        return ZTOP, VP, VS, RH
    #------------------
    def startmodel(self, sigH = 0.01, sigVS = 0.1, smoothVS = 0., lockbot = False):
        mapr  = self.__call__(self.ZTOP, self.VP, self.VS, self.RH)
        H = self.ZTOP[1:] - self.ZTOP[:-1]
        siglogH = sigH / abs(H)
        siglogVS = sigVS / self.VS

        if lockbot:
            #cannot lock sigH   
            assert len(siglogVS) > 1
            siglogH[-1] = siglogVS[-1] = 1.e-20


        CMinv = np.concatenate((siglogH, siglogVS)) ** -2.

        if smoothVS:
            N = len(self.VS)
            rho = np.eye(len(mapr))
            r    = np.arange(N)
            I, J = np.meshgrid(r, r)
            D = abs(I - J)
            rho[N-1:2*N-1, N-1:2*N-1] = np.exp(-D / smoothVS) #smoothing distance in cell number
            if False:
                plt.figure()
                gca().matrixview(rho)
                showme()
            VX, VY = np.meshgrid(1. / CMinv, 1. / CMinv)
            CM = rho * np.sqrt(VX * VY)
            CMinv = np.linalg.inv(CM)

        return mapr, CMinv
#------------------
class Parameterizer_logVS_PRlaw_RHlaw(Parameterizer_logH_logVS_PRlaw_RHlaw):
    """logH and logVS are parameters,
       VP is recovered as a function of VS using VP/VS = f(depth)
       RH is recovered as a function of VP using RHO = f(VP)
    """
    def __call__(self, ZTOP, VP, VS, RH):
        "ignore VP RH"
        return np.log(VS)
    #------------------
    def inv(self, m):
        logVS  = m

        VS   = np.exp(logVS)
        ZTOP = self.ZTOP
        ZMID = np.concatenate((0.5 * (ZTOP[1:] + ZTOP[:-1]), [1.1 * ZTOP[-1]]))

        VP = VS * self.prlaw(ZMID)
        RH = self.rhlaw(VP)

        return ZTOP, VP, VS, RH
    #------------------
    def startmodel(self, sigVS = 0.1, smoothVS = 0., lockbot = False):
        mapr  = self.__call__(self.ZTOP, self.VP, self.VS, self.RH)

        siglogVS = sigVS / self.VS

        if lockbot:
            #cannot lock sigH   
            assert len(siglogVS) > 1
            siglogVS[-1] = 1.e-20

        CMinv = siglogVS ** -2.

        if smoothVS:
            N = len(self.VS)
            rho = np.eye(len(mapr))
            r    = np.arange(N)
            I, J = np.meshgrid(r, r)
            D = abs(I - J)
            rho = np.exp(-D / smoothVS) #smoothing distance in cell number
            if False:
                plt.figure()
                gca().matrixview(rho)
                showme()
            VX, VY = np.meshgrid(1. / CMinv, 1. / CMinv)
            CM = rho * np.sqrt(VX * VY)
            CMinv = np.linalg.inv(CM)

        return mapr, CMinv
#------------------
class Parameterizer_logH_PRlaw_RHlaw(Parameterizer_logH_logVS_PRlaw_RHlaw):
    """logH and logVS are parameters,
       VP is recovered as a function of VS using VP/VS = f(depth)
       RH is recovered as a function of VP using RHO = f(VP)
    """

    #------------------
    def __call__(self, ZTOP, VP, VS, RH):
        "ignore VP RH"
        H  = ZTOP[1:] - ZTOP[:-1]
        return np.log(H)
    #------------------
    def inv(self, m):
        logH  = m

        H    = np.exp(logH)
        ZTOP = np.concatenate(([0.], H.cumsum()))
        ZMID = np.concatenate((0.5 * (ZTOP[1:] + ZTOP[:-1]), [1.1 * ZTOP[-1]]))

        VS = self.VS
        VP = VS * self.prlaw(ZMID)
        RH = self.rhlaw(VP)

        return ZTOP, VP, VS, RH
    #------------------
    def startmodel(self, sigH = 0.1, lockbot = False):
        mapr  = self.__call__(self.ZTOP, self.VP, self.VS, self.RH)
        H = self.ZTOP[1:] - self.ZTOP[:-1]
        siglogH = sigH / abs(H)

        if lockbot:
            raise Exception('')

        CMinv = siglogH ** -2.
        return mapr, CMinv




###################
def displininv(s96, mod96, \
        sigH       = 0.1, 
        sigVS      = 0.1, 
        sigPR      = 0.1, 
        sigRH      = 0.1, 
        smoothVS   = 0., 
        smoothPR   = 0., 
        smoothRH   = 0., 
        damping    = 1.0, 
        imax       = 50, 
        interruptconditon = 0.01, 
        lockbot    = False):
    """
    s96    = target data, filename at SURF96 format or string
    mod96  = starting mode, filename at MOD96 format or string
    sigH   = uncertainty on prior layer thickness in km
    sigVS  = uncertainty on VS in km/s
    ...
    """
    #--------------
    D = makedatacoder(s96, which = Datacoder_log)
    dobs, CDinv = D.target()

    #--------------
    if sigH and sigVS and sigPR and sigRH:
        #print "  invert all"
        P = makeparameterizer(mod96, 
            which = Parameterizer_logH_logVS_logPR_logRH)
        mapr, CMinv = P.startmodel(sigH = sigH, sigVS = sigVS, sigPR = sigPR, sigRH = sigRH, smoothVS = smoothVS, smoothPR = smoothPR, smoothRH = smoothRH, lockbot = lockbot)

    elif sigH and sigVS and not sigPR and not sigRH:
        #print "  lock PR and RH"
        P = makeparameterizer(mod96, 
            which = Parameterizer_logH_logVS)
        mapr, CMinv = P.startmodel(sigH = sigH, sigVS = sigVS, smoothVS = smoothVS, lockbot = lockbot)

    elif not sigH and sigVS and not sigPR and not sigRH:
        #print "  lock H, PR and RH"
        P = makeparameterizer(mod96, 
            which = Parameterizer_logVS)
        mapr, CMinv = P.startmodel(sigVS = sigVS, smoothVS = smoothVS, lockbot = lockbot)

    elif sigH and not sigVS and not sigPR and not sigRH:
        #print "  lock VS, PR and RH"
        P = makeparameterizer(mod96, 
            which = Parameterizer_logH)
        mapr, CMinv = P.startmodel(sigH = sigH, lockbot = lockbot)

    elif sigH and sigVS and sigPR and not sigRH:
        raise NotImplementedError('lock RH')
    elif not sigH and sigVS and sigPR and not sigRH:
        raise NotImplementedError("lock H and RH")
    else:
        raise NotImplementedError('')

    #--------------
    g = Theory(P, D)

    #--------------
    L = LinInv(\
        g         = g, 
        d         = dobs, 
        CDinv     = CDinv,
        mapr      = mapr, 
        CMinv     = CMinv,
        damping   = damping,
        eps       = 1.e-2, 
        method    = 1)
    
    #--------------
    mt, dt, chi2red = L.upd()
    CHIs = [chi2red]
    ZTOP, VP, VS, RH = P.inv(mt)
    values = D.inv(dt)
    yield chi2red, (ZTOP, VP, VS, RH), (D.waves, D.types, D.modes, D.freqs, values)

    #--------------
    for n in xrange(imax):
        try:
            mt, dt, G = L.grad()
            mb, db = L.inv()
            mt, dt, chi2red = L.upd()
        except KeyboardInterrupt: raise
        except Exception as e: 
            print ("\n", e)
            #raise
            break

        CHIs.append(chi2red)
        ZTOP, VP, VS, RH = P.inv(mt)
        values = D.inv(dt)
        yield chi2red, (ZTOP, VP, VS, RH), (D.waves, D.types, D.modes, D.freqs, values)

        if len(CHIs) > 5 and (np.abs(CHIs[-5:] / CHIs[-1] - 1.) <= interruptconditon).all():
            break
###################
def display(displininvout):
    """
    input = list(displininv(*args, **kwargs))
    """
    rd = resultdisplay()
    #target = displininvout.pop(0) #first item is the target data
    #rd.plotdisp(color =  "k", linewidth = 3, alpha = 1, *target)

    N = len(displininvout)
    for n, (chi2red, model, data) in enumerate(displininvout):
        clr = value2color(n / float(N), cmap = tej())
        if not n:
            rd.plotmodel(color = clr, linewidth = 3, alpha = 1, *model)
            rd.plotdisp(color =  clr, linewidth = 3, alpha = 1, *data)
        else:
            rd.plotmodel(color = clr, linewidth = 1, alpha = 1.0, *model)
            rd.plotdisp(color  = clr, linewidth = 1, alpha = 1.0, *data)
        rd.axconv.plot(n, chi2red, "o", color=  clr)
    rd.plotmodel(color = clr, linewidth = 3, alpha = 1.0, *model)
    rd.plotdisp(color  = clr, linewidth = 3, alpha = 1.0, *data)

    rd.grid()
    rd.tick()
    return rd
###################

###################
if __name__ == "__main__":

    o = list(displininv(\
        s96      = "/home/max/progdat/CPiS/EarthModel/bidon0.surf96",
        mod96    = "/home/max/progdat/CPiS/EarthModel/Soultz.rho.mod", 
        sigH       = 0.1, 
        sigVS      = 0.1, 
        sigPR      = 0.1, 
        sigRH      = 0.1, 
        smoothVS   = 0., 
        smoothPR   = 0., 
        smoothRH   = 0., 
        damping    = 1.0, 
        imax       = 50, 
        interruptconditon = 0.01, 
        lockbot    = False))
    rd = display(o)
    rd.plotm96("/home/max/progdat/CPiS/EarthModel/bidon0.mod96", color = 'k', linewidth = 3, alpha = 1.0)
    rd.plots96("/home/max/progdat/CPiS/EarthModel/bidon0.surf96", color = 'k', linewidth = 3, alpha = 1.0)
    
    showme()


#def displininv(s96, mod96, \
#        sigH       = 0.1, 
#        sigVS      = 0.1, 
#        smoothVS   = 0., 
#        damping    = 1.0, 
#        imax       = 50, 
#        interruptconditon = 0.01, 
#        prlaw      = "1.0335 * np.exp(-z / 0.5408) + 1.7310",
#        rhlaw      = "1.74 * vp ** 0.25"):
#    """
#    s96    = target data, filename at SURF96 format or string
#    mod96  = starting mode, filename at MOD96 format or string
#    sigH   = uncertainty on prior layer thickness in km
#    sigVS  = uncertainty on VS in km/s
#    ...
#    """
#    if "grt1" in prlaw.lower(): prlaw = eval('lambda z : 1.0335 * np.exp(-z / 0.5408) + 1.7310')
#    else:                       prlaw = eval('lambda z : %s' % prlaw)
#    if "gardner" in rhlaw.lower(): rhlaw = eval('lambda vp : 1.74 * vp ** 0.25')
#    else:                          rhlaw = eval('lambda vp : %s' % rhlaw)
#    #--------------
#    D = makedatacoder(s96, which = Datacoder_log)
#    dobs, CDinv = D.target()
#    #yield (D.waves, D.types, D.modes, D.freqs, D.values, D.dParameterizer_logH_logMU_logK_logRHvalues)

#    #--------------
#    if sigH and sigVS:
#        P = makeparameterizer(mod96, 
#            which = Parameterizer_logH_logVS_PRlaw_RHlaw, 
#            prlaw = prlaw, 
#            rhlaw = rhlaw)
#        mapr, CMinv = P.startmodel(sigH = sigH, sigVS = sigVS, smoothVS = smoothVS)
#    elif sigH == 0. and sigVS:
#        P = makeparameterizer(mod96, 
#            which = Parameterizer_logVS_PRlaw_RHlaw, 
#            prlaw = prlaw, 
#            rhlaw = rhlaw)
#        mapr, CMinv = P.startmodel(sigVS = sigVS, smoothVS = smoothVS)
#    elif sigVS == 0. and sigH:
#        P = makeparameterizer(mod96, 
#            which = Parameterizer_logH_PRlaw_RHlaw, 
#            prlaw = prlaw, 
#            rhlaw = rhlaw)
#        mapr, CMinv = P.startmodel(sigH = sigH)
#    else: raise Exception('')
#    
#    #--------------
#    g = Theory(P, D)

#    #--------------
#    L = LinInv(\
#        g         = g, 
#        d         = dobs, 
#        CDinv     = CDinv,
#        mapr      = mapr, 
#        CMinv     = CMinv,
#        damping   = damping,
#        eps       = 1.e-3, 
#        method    = 1)
#    
#    #--------------
#    mt, dt, chi2red = L.upd()
#    CHIs = [chi2red]
#    ZTOP, VP, VS, RH = P.inv(mt)
#    values = D.inv(dt)
#    yield chi2red, (ZTOP, VP, VS, RH), (D.waves, D.types, D.modes, D.freqs, values)

#    #--------------
#    for n in xrange(imax):
#        try:
#            mt, dt, G = L.grad()
#            mb, db = L.inv()
#            mt, dt, chi2red = L.upd()
#        except KeyboardInterrupt: raise
#        except Exception as e: 
#            #print "\n", e
#            #raise
#            break

#        CHIs.append(chi2red)
#        ZTOP, VP, VS, RH = P.inv(mt)
#        values = D.inv(dt)
#        yield chi2red, (ZTOP, VP, VS, RH), (D.waves, D.types, D.modes, D.freqs, values)

#        if len(CHIs) > 5 and (np.abs(CHIs[-5:] / CHIs[-1] - 1.) <= interruptconditon).all():
#            break

#    mt, dt, chi2red = L.upd()
#    rd.plotmodel(color = "r", linewidth = 3, alpha = 1, *P.inv(mt))
#    rd.plotdisp(D.waves, D.types, D.modes, D.freqs, D.inv(dt), None, color = "r", alpha = 1.0, linewidth = 3)

#    print init
#    print target
#    for m in Ms:
#        print m
#    for d in Ds:
#        print d
#    for chi2red in CHIs:
#        print chi2red



#    from labex import timer
#    truemodel = depthmodel_from_mod96("/home/max/progdat/CPiS/EarthModel/bidon0.mod96")
#    D = makedatacoder("/home/max/progdat/CPiS/EarthModel/bidon0.surf96")

#    if False:
#        P = makeparameterizer("/home/max/progdat/CPiS/EarthModel/Soultz.rho.mod", 
#                which = Parameterizer_logH_logVS_logPR_logRH)
#        mapr, CMinv = P.startmodel(sigH = 0.01, sigVS = 0.3, sigPR = 0.2, sigRH = 0.1, smoothVS = 0., smoothPR = 0., smoothRH = 0.)
#    elif False:
#        P = Parameterizer_logH_logVS(\
#                ZTOP = np.linspace(0., 3., 5),
#                VP   = 2.0 * np.ones(5),        
#                VS   = 1.0 * np.ones(5),
#                RH   = 2.3 * np.ones(5))
#        mapr, CMinv = P.startmodel(sigH = 0.1, sigVS = 0.3, smoothVS = 0.)
#    elif False:
#        P = Parameterizer_logVS(\
#                ZTOP = np.linspace(0., 3., 5),
#                VP   = 2.0 * np.ones(5),        
#                VS   = 1.0 * np.ones(5),
#                RH   = 2.3 * np.ones(5))
#        mapr, CMinv = P.startmodel(sigVS = 0.3, smoothVS = 10.)
#    elif False:
#        P = makeparameterizer("/home/max/progdat/CPiS/EarthModel/Soultz.rho.mod", 
#                which = Parameterizer_logH)
#        mapr, CMinv = P.startmodel(sigH = 0.1)
#    elif False:
#        P = Parameterizer_logH_logVS_PRlaw_RHlaw(\
#                ZTOP = np.linspace(0., 3., 4),
#                VP   = None,        
#                VS   = 1.0 * np.ones(4),
#                RH   = None)
#        mapr, CMinv = P.startmodel(sigH = 0.1, sigVS = 0.1, smoothVS = 2.)
#    elif False:
#        P = Parameterizer_logVS_PRlaw_RHlaw(\
#                ZTOP = np.linspace(0., 3., 4),
#                VP   = None,        
#                VS   = 1.0 * np.ones(4),
#                RH   = None)
#        mapr, CMinv = P.startmodel(sigVS = 0.1, smoothVS = 2.)
#    elif True:
#        P = makeparameterizer("/home/max/progdat/CPiS/EarthModel/Soultz.rho.mod", 
#                which = Parameterizer_logH_PRlaw_RHlaw)
#        mapr, CMinv = P.startmodel(sigH = 0.1)

#    dobs, CDinv = D.target()
#    g = Theory(P, D)

#    L = LinInv(\
#        g         = g, 
#        d         = dobs, 
#        CDinv     = CDinv,
#        mapr      = mapr, 
#        CMinv     = CMinv,
#        damping   = 2.,
#        eps       = 1.e-3, 
#        method    = 1)
#    
#    rd = resultdisplay()
#    rd.plotdisp(D.waves, D.types, D.modes, D.freqs, D.values, D.dvalues, color = "k", alpha = 1, linewidth = 3)

#    mt, dt, chi2red = L.upd()
#    rd.plotmodel(color = "r", linewidth = 3, alpha = 1, *P.inv(mt))
#    rd.plotdisp(D.waves, D.types, D.modes, D.freqs, D.inv(dt), None, color = "r", alpha = 1.0, linewidth = 3)

#    for n in xrange(10):
#        with timer(str(n)):
#            mt, dt, G = L.grad()
#        with timer(str(n)):
#            mb, db = L.inv()
#        with timer(str(n)):
#            mt, dt, chi2red = L.upd()
#        rd.plotmodel(color = "b", linewidth = 1, alpha = 0.2, *P.inv(mt))
#        rd.plotdisp(D.waves, D.types, D.modes, D.freqs, D.inv(dt), None, color = "b", alpha = 0.2, linewidth = 1)
#        print mt
#    rd.plotmodel(color = "g", linewidth = 3, alpha = 1, *P.inv(mt))
#    rd.plotdisp(D.waves, D.types, D.modes, D.freqs, D.inv(dt), None, color = "g", alpha = 1, linewidth = 3)
#    rd.tick()
#    rd.grid()

#    rd.plotmodel(truemodel.vp.z,
#                 truemodel.vp.values, 
#                 truemodel.vs.values, 
#                 truemodel.rh.values, 
#                 color = "k", linewidth = 3, alpha = 1)

#    showme()

##------------------
#class G(object):
#    h   = 0.005
#    dcl = 0.005
#    dcr = 0.005
#    #------------------
#    def __init__(self, waves, types, modes, freqs, values, dvalues):
#        assert np.all(~np.isnan(values))
#        assert np.all([len(w) == len(waves) for w in types, modes, freqs, values, dvalues])
#        self.waves = waves     #parameters for forward problem
#        self.types = types     #parameters for forward problem
#        self.modes = modes     #parameters for forward problem
#        self.freqs = freqs     #parameters for forward problem
#        self.values = values   #target data
#        self.dvalues = dvalues #target data
#    #------------------
#    def __call__(self, m):
#        ZTOP, VP, VS, RH = logH_logMU_logK_logRH(m)
#        values = srfdis17(ZTOP, VP, VS, RH, \
#                self.waves, self.types, self.modes, self.freqs,
#                self.h, self.dcl, self.dcr)
#        return np.log(values)
##------------------
#class G_from_s96(G):
#    def __init__(self, s96):
#        s = surf96reader(s96)
#        waves, types, modes, freqs, values, dvalues = s.wtmfvd()
#        G.__init__(self, waves, types, modes, freqs, values, dvalues)
##------------------
#class Mapr_from_mod96(mod96):
#    dm = depthmodel_from_mod96(mod96)
#    ztop = dm.vp.z
#    vp   = dm.vp.value
#    vs   = dm.vs.value
#    rh   = dm.rh.value
#    mu   = rh * vs ** 2.
#    K    = rh * vp ** 2. - 4. * mu / 3.
#    logH  = np.log(ztop[1:] - ztop[:-1])
#    logMU = np.log(mu)
#    logK  = np.log(K)
#    logRH = np.log(rh)
#    return np.concatenate((logH, logMU, logK, logRH))
##------------------
#if __name__ == "__main__":
#    g = G_from_s96("datatest.surf96")
#    print g.waves
#    print g.types
#    print g.modes
#    print g.freqs
#    print g.values



#    #------------------
#    class Parameterizer_logHMUK(Parameterizer_logH_logMU_logK_logRH):
#        """
#        """
#        def __call__(self, ZTOP, VP, VS, RH):
#            "ingore input RH"
#            assert len(VP) == len(VP) == len(ZTOP) == len(self.RH)
#            H  = ZTOP[1:] - ZTOP[:-1]
#            MU = self.RH * VS ** 2.
#            K  = self.RH * VP ** 2. - 4. * MU / 3.
#            return np.log(np.concatenate((H, MU, K)))
#        #------------------
#        def inv(self, m):
#            assert not (len(m) + 1) % 3
#            nlayer = (len(m) + 1) / 3
#            logH   = m[              : 1 * nlayer -1]
#            logMU  = m[1 * nlayer - 1: 2 * nlayer -1]
#            logK   = m[2 * nlayer - 1: 3 * nlayer -1]

#            H    = np.exp(logH)
#            MU   = np.exp(logMU)
#            K    = np.exp(logK)

#            ZTOP = np.concatenate(([0.], H.cumsum()))
#            VP   = np.sqrt((K + 4. * MU / 3.) / self.RH)
#            VS   = np.sqrt(MU / self.RH)
#            return ZTOP, VP, VS, self.RH
#        #------------------
#        def startmodel(self, sigz = 0.1, sigvs = 1.0, sigpr = 0.1):
#            mapr  = self.__call__(self.ZTOP, self.VP, self.VS, self.RH)
#            H = self.ZTOP[1:] - self.ZTOP[:-1]

#            siglogH  = 2. * sigz / abs(H)
#            siglogMU = 2. * sigvs / self.VS
#            PR = self.VP / self.VS
#            siglogK = siglogMU + abs(2. * PR / (PR ** 2. - 4. / 3.)) * sigpr
#            CMinv = np.concatenate((siglogH, siglogMU, siglogK)) ** -2.
#            return mapr, CMinv
#    #------------------
#    class Parameterizer_logH_logVS(Parameterizer_logH_logMU_logK_logRH):
#        def __call__(self, ZTOP, VP, VS, RH):
#            "ingore input VP and RH"
#            assert len(VS) == len(ZTOP) == len(self.RH)
#            H  = ZTOP[1:] - ZTOP[:-1]
#            return np.log(np.concatenate((H, VS)))
#        #------------------
#        def inv(self, m):
#            assert not (len(m) + 1) % 2
#            nlayer = (len(m) + 1) / 2
#            logH   = m[              : 1 * nlayer -1]
#            logVS  = m[1 * nlayer - 1: 2 * nlayer -1]

#            H    = np.exp(logH)
#            VS   = np.exp(logVS)

#            ZTOP = np.concatenate(([0.], H.cumsum()))
#            VP   = VS * (self.VP / self.VS)

#            return ZTOP, VP, VS, self.RH
#        #------------------
#        def startmodel(self, sigz = 0.1, sigvs = 1.0, smoothvs = 10.):
#            mapr  = self.__call__(self.ZTOP, self.VP, self.VS, self.RH)
#            H = self.ZTOP[1:] - self.ZTOP[:-1]
#            siglogH  = 2. * sigz / abs(H)
#            siglogVS = sigvs / abs(self.VS)
#            CMinv = np.concatenate((siglogH, siglogVS)) ** -2.

#            #####
#            if smoothvs:
#                rho = np.eye(len(mapr))
#                r    = np.arange(len(self.VS))
#                I, J = np.meshgrid(r, r)
#                D = abs(I - J)
#                rho[len(H):, len(H):] = np.exp(-D / smoothvs) #smoothing distance in cell number
#                VX, VY = np.meshgrid(1. / CMinv, 1. / CMinv)
#                CM = rho * np.sqrt(VX * VY)
#                CMinv = np.linalg.inv(CM)
#            #####
#            return mapr, CMinv
