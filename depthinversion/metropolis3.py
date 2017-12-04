import numpy as np
from numpy import log
import time

#########################################################
class LogUni(object):
    """logarithm of a 1D, uniform, non normalized pdf (pseudo-uniform)
        I nead a pdf approaching a uniform distribution
        I want "bad" models to keep comparable

        f(x) =  | 0                          if     vinf < x < vsup
                | -1/2 * (x - vinf) / vstd   if     x < vinf
                | -1/2 * (x - vsup) / vstd   if     x > vsup
                | nanpenalty  or exception   if     x is nan or inf
                where vstd is (vsup - vinf) / k, k > 1
    """
    def __init__(self, vinf, vsup, k=1000.0, nanbehavior=0):
        """

        :param vinf: float, lowest bound for uni pdf
        :param vsup: float, highest bound for uni pdf
        :param k: float, sharpness coefficients for the edges of the non-zero probability area
        :param nanbehavior: int,
            0 means raise an exception in case of nans
            1 means add a gentle penalty to models which include nans
            2 means add a severe penalty to models which include nans
        """
        #vinf, vsup = floats
        assert not np.any(np.isnan([vinf, vsup, k]))
        assert not np.any(np.isinf([vinf, vsup, k]))
        assert vinf < vsup
        assert k > 1.
        self.vinf, self.vsup = vinf, vsup
        self.vmean = 0.5 * (self.vinf + self.vsup)
        self.std = (self.vsup - self.vinf) / float(k)
        self.nanbehavior = nanbehavior
        if self.nanbehavior == 0:
            self.raiseifnan  = True

        elif self.nanbehavior == 1:
            # case 1: return a constant penalty, for ND pdfs, the penalty will be proportionnal to the number of nans
            self.raiseifnan = False
            self.nanpenalty = self.call1(self.vinf - 0.1 * (self.vsup - self.vinf)) #take the pdf value away from the 0 area

        elif self.nanbehavior == 2:
            # case 2: return severe penalty
            self.raiseifnan = False
            self.nanpenalty = self.call1(self.vinf - 10000. * (self.vsup - self.vinf))  # take the pdf value away from the 0 area
        else: raise ValueError('no such nan behavior implemented')

    def calln(self, v):
        """v is an array with unefined length"""
        y = np.zeros_like(v)
        #-----------------------
        GOOD = ~(np.isnan(v) | np.isinf(v))
        if not GOOD.all():
            if self.raiseifnan: raise Exception('pdf got innapropriate values')
            else: y[~GOOD] = self.nanpenalty


        #-----------------------
        I = J = np.zeros_like(GOOD)

        I[GOOD] = (v[GOOD] < self.vinf)
        y[I] += -0.5 * ((v[I] - self.vinf) / self.std) ** 2. 

        J[GOOD] = (self.vsup < v[GOOD])
        y[J] += -0.5 * ((v[J] - self.vsup) / self.std) ** 2. 
        #-----------------------
        return y

    def call1(self, v):
        """v is a scalar"""
        assert not hasattr(v, "__len__")
        #return self.calln(np.asarray([v]))[0]  #too slow
        if np.isnan(v) or np.isinf(v):
            if self.raiseifnan: raise Exception('pdf got innapropriate value')
            else: return self.nanpenalty

        if v < self.vinf:
            return -0.5 * ((v - self.vinf) / self.std) ** 2.
        elif self.vsup < v:
            return -0.5 * ((v - self.vsup) / self.std) ** 2.
        return 0.

    def __call__(self, v): return self.call1(v)
#########################################################
class LogGauss(object):
    """1D truncated gaussian, take the product of a gaussian and a pseudo-uniform distribution
        the nan or inf behavior is handled by the pseudo uniform pdf

        f(x, args) = LogUni(x, args)  -1/2 * ((x - vmean) / vstd) ** 2.
    """
    def __init__(self, vmean, vunc, vinf, vsup, k = 1000., nanbehavior = 0):
        assert not np.any(np.isnan([vmean, vunc, vinf, vsup, k]))
        assert not np.any(np.isinf([vmean, vunc, vinf, vsup, k]))
        assert vinf <= vmean <= vsup
        assert k > 0.
        self.vinf, self.vsup = vinf, vsup
        self.vmean, self.vunc = vmean, vunc 
        self.luni = LogUni(vinf, vsup, k, nanbehavior)

    def calln(self, v):
        """v is an array with unefined length"""
        y = self.luni.calln(v) #handles nans or infs according to nanbehavior

        K = ~(np.isnan(v) | np.isinf(v))
        y[K] += -0.5 * ((v[K] - self.vmean) / self.vunc) ** 2. 

        return y

    def call1(self, v):
        """v is a scalar"""
        y = self.luni.call1(v)
        if np.isnan(v) or np.isinf(v): return y
        return y + -0.5 * ((v - self.vmean) / self.vunc) ** 2.

    def __call__(self, v): return self.call1(v)
#########################################################
class LogUniND(object):
    """N uniform laws, assume independent variables, sum of LogUni objects"""
    def __init__(self, vinfs, vsups, k = 1000., nanbehavior = 0):
        vinfs, vsups = [np.asarray(_, float) for _ in [vinfs, vsups]]
        assert len(vinfs) == len(vsups)
        self.N = len(vinfs)
        self.ls = [LogUni(vinf, vsup, k, nanbehavior) for vinf, vsup in zip(vinfs, vsups)]

    def callpoints(self, points):
        #points = 2D : raws = list of points, columns = dimensions
        y = np.zeros_like(points.shape[0], float)
        for n, l in enumerate(self.ls):
            y += l(points[:, n])
        return y

    def callargs(self, *args):
        #args like (X0, X1, X2, ..., XM-1) where Xi are arrays with same shapes
        assert len(args) == self.N
        for arg in args: assert arg.shape == args[0].shape

        y = np.zeros_like(arg)
        for n, (l, arg) in enumerate(zip(self.ls, args)):
            y += l.calln(arg)
        return y

    def call1(self, model):
        #model like [x0, x1, x2, ..., xM-1] where xi are scalars
        assert len(model) == self.N
        y = sum([l(model[j]) for j, l in enumerate(self.ls)])
        return y

    def __call__(self, model):
        #model like [x0, x1, x2, ..., xM-1] where xi are scalars
        return self.call1(model)
#########################################################
class LogGaussND(LogUniND):
    """ND truncated gaussians, no covariance, sum of LogGauss objects"""
    def __init__(self, vmeans, vuncs, vinfs, vsups, k = 1000., nanbehavior = 0):
        vmeans, vuncs, vinfs, vsups = [np.asarray(_, float) for _ in [vmeans, vuncs, vinfs, vsups]]
        assert len(vinfs) == len(vsups) == len(vmeans) == len(vuncs)
        self.N = len(vinfs)
        self.ls = [LogGauss(vmean, vunc, vinf, vsup, k, nanbehavior) \
                        for vmean, vunc, vinf, vsup in zip(vmeans, vuncs, vinfs, vsups)]
#########################################################
class LogGaussNDCov(LogUniND):
    """ND truncated gaussians, with covariance"""
    def __init__(self, vmeans, vuncs, vinfs, vsups, rho, k = 1000., nanbehavior = 0):
        vmeans, vuncs, vinfs, vsups, rho = [np.asarray(_, float) for _ in [vmeans, vuncs, vinfs, vsups, rho]]
        assert len(vinfs) == len(vsups) == len(vmeans) == len(vuncs)
        assert rho.shape == (len(vinfs), len(vinfs))
        assert np.all((-1.0 <= rho) & (rho <= 1.0))
        self.N = len(vinfs)

        self.luni = LogUniND(vinfs, vsups, k, nanbehavior)
        sX, sY = np.meshgrid(vuncs, vuncs)
        self.CDinv = np.linalg.inv(rho * sX * sY)
        self.vmeans = np.asarray(vmeans)

    def callargs(self, *args):
        #args like (X0, X1, X2, ..., XM-1) where Xi are arrays with same shapes
        y = self.luni.callargs(*args) #includes penalty for nans or infs
        for n, model in enumerate(zip(*[arg.flat[:] for arg in args])):
            y.flat[n] += self.call1(model)
        return y
    
    def call1(self, model):
        #model like [x0, x1, x2, ..., xM-1] where xi are scalars
        M = np.asarray(model) - self.vmeans
        y = -.5 * np.dot(np.dot(M, self.CDinv), M.T) + self.luni.call1(model)
        return y

    def __call__(self, model):
        #model like [x0, x1, x2, ..., xM-1] where xi are scalars
        return self.call1(model)
#########################################################
def metropolis(M0, MSTD, G, ND, logRHOD, logRHOM, nkeep, normallaw=np.random.randn, unilaw=np.random.rand,
               chainid=1, HL=100, IK0=0.25, MPMIN=1.e-6, MPMAX=1e6, adjustspeed=0.01, nofail=False, debug=False):
    """
    input :
        M0      = starting model, np.ndarray
        MSTD    = markov std array, np.ndarray
        G       = theory function, callable : np.ndarray (model) -> np.ndarray (data) 
        ND      = integer, number of dimensions in the data space
        logRHOD = data prior pdf, callable  : np.ndarray (data) -> float (log(pdf(data)))
        logRHOM = model prior pdf, callable : np.ndarray (model) -> float (log(pdf(model)))
        nkeep   = number of models to keep
        HL      = history length, used to adjust the proposal according to the recent acceptance ratio and recent std in each dimension
        IK0     = target acceptance ratio
        adjustspeed = a coefficient controlling the update of the masterproposal term according to the disntance between
                      the current acceptance ratio and the target one (IK0), 0 means no adjustment, use < 1.
        nofail  = if all call to G are failing, I will suspect a programming error in the theory, rather than unsuited models.
                  use nofail = True to avoid raising an exception
        debug   = if True, raise in case of theory failure, and display the error message
    output : 
        models  = 2D array, each line is a model retained, each model appears only once
        datas   = 2D array, corresponding data arrays, nans means theory failure
        weights = 1D array, number of times each model was retained
        llks    = 1D array, log likelyhood associated to each model

        !! use !!
            models.repeat(weights, axis = 0) if you need each model to be repeated the correct number of times (also for datas, llks)
        !! notes !!
            * for parallel computing avoid passing functions G, logRHOD and logRHOM through pipes for pickling reasons
            * use thread safe random functions, numpy onces are not thread safe
            * failure during call to G are treated as non fatal, if the theory fails to produce data, the likelyhood is penalyzed so that the chain should run 
              away from that zones. If all models are failing, the chain will sample the prior pdf only (with very low llks)                              
    #-------------------
    variables:
        ntest    = number of tested models since chain start
        nfail    = number of models for which data could not be computed since chain start (the chain still samples the prior pdf)
        nkept    = number of distinct(!) models retained since chain start
        nstay    = number of subsequent rejection (= current weight of the last generated model)
        icurrent = the position in array models of the last generated model

    """
    nfail = 0
    #----      
    def call(M, nfail):
        try:
            D = G(M)
            L = logRHOM(M) + logRHOD(D)
        except KeyboardInterrupt: raise
        except Exception as e:
            if debug:
                raise
                #helps user to understand why the theory fails, still not raise an exception
                print e
            #could not compute data array, still sample the prior pdf to converge toward successfull regions, add constant penalty
            nfail += 1
            D = np.zeros(ND, float) * np.nan  
            L = logRHOM(M) + -1e20 #bad models remain comparable, try to guide the chain away from these zones based on prior info
        return D, L, nfail
    #----        
    MI = M0
    DI, LI, nfail = call(MI, nfail)
    #----
    ntest = nkept = nstay = 0
    MP    = 1.0 #master proposal, used to scale the markov chain proposal stds
    Ikept = np.zeros(HL, bool) #chain status history
    start = time.time()
    lasttime, lastntest = start, 0


    #---- chain data
    #store every distinct model and data of the chain 
    models  = np.empty((nkeep+1, len(MSTD)), float)
    datas   = np.empty((nkeep+1, len(DI)), float)
    llks    = np.zeros(nkeep+1, float) + -np.inf
    weights = np.zeros(nkeep+1, int)

    icurrent            = 0
    models[icurrent, :] = MI
    datas[icurrent, :]  = DI
    llks[icurrent]      = LI
    weights[icurrent]  += 1
    #---- 
    while nkept < nkeep:
        ntest += 1
        #----------------------
        if not ntest % HL and ntest:
            #reevaluate stats and master proposal
            IK  = Ikept.sum() / float(HL)  #instantaneous keep ratio, gives the recent acceptance of the chain
            AK  = nkept / float(ntest)     #average keep ratio
            MP *= 1. - adjustspeed * (IK0 - IK)   #master proposal, used to scale the markov chain proposal stds
            MP = np.clip(MP, MPMIN, MPMAX)
            t = time.time()
            AS  = ntest / (t - start)                  #Average test speed
            IS  = (ntest - lastntest) / (t - lasttime) #Instantaneous test speed
            lasttime, lastntest = t, ntest

        #----------------------
        if not ntest % (10 * HL) and ntest:
            #reevaluate proposal stds
            I = np.arange(np.max([0, icurrent - 10 * HL]), icurrent + 1)
            MSTD = np.std(models[I, :].repeat(weights[I], axis = 0), axis = 0)
        #----------------------
        if nfail >= nkeep and nfail >= ntest and not nofail:
            #all attempts to call G failed, it might be due to a programming error...
            print ("chainid %5d ntest %5d nfail %5d nkept %5d nstay %5d IK %5.2f AK %5.2f MP %.2f AS %6.2f/s IS %6.2f/s LI %f, PRESUMED ERROR IN THEORY" % (chainid, ntest, nfail, nkept, nstay, IK, AK, MP, AS, IS, LI))
            while icurrent:
                G(models[icurrent, :]) #run G to reproduce the error message
                icurrent -= 1
            raise Exception('should never occur')
        #----------------------
        if nstay >= nkeep:
            #stuck signal
            print ("chainid %5d ntest %5d nfail %5d nkept %5d nstay %5d IK %5.2f AK %5.2f MP %.2f AS %6.2f/s IS %6.2f/s LI %f, STUCK" % (chainid, ntest, nfail, nkept, nstay, IK, AK, MP, AS, IS, LI))
            return models[:icurrent, :], datas[:icurrent, :], weights[:icurrent], llks[:icurrent]
        #----------------------
        if not ntest % (10 * HL) and ntest: #True
            #report stats
            print ("chainid %5d ntest %5d nfail %5d nkept %5d nstay %5d IK %5.2f AK %5.2f MP %.2f AS %6.2f/s IS %6.2f/s LI %f" % (chainid, ntest, nfail, nkept, nstay, IK, AK, MP, AS, IS, LI))
            #print (MP, " ".join(['%.4f' % _ for _ in MSTD / MSTD.max()]))
        #----------------------
        MII = MI + MP * MSTD * normallaw(len(M0))
        DII, LII, nfail = call(MII, nfail)
        #----------------------
        if LII > LI:
            MI, DI, LI = MII, DII, LII
            nkept += 1
            nstay  = 0
            Ikept[ntest % HL] = True

            icurrent += 1
            models[icurrent, :] = MI
            datas[icurrent, :] = DI
            llks[icurrent] = LI
            weights[icurrent] += 1
            #yield False, MI, DI, LI #improvement
        #----------------------
        else:
            r =  unilaw()
            if r == 0. or log(r) < (LII - LI):
                MI, DI, LI = MII, DII, LII
                nkept += 1
                nstay  = 0
                Ikept[ntest % HL] = True

                icurrent += 1
                models[icurrent, :] = MI
                datas[icurrent, :] = DI
                llks[icurrent] = LI
                weights[icurrent] += 1

                #yield False, MI, DI, LI #decay
            else:
                nstay += 1
                Ikept[ntest % HL] = False
                weights[icurrent] += 1
                #yield True, MI, DI, LI #stay
    print ("chainid %5d ntest %5d nfail %5d nkept %5d nstay %5d IK %5.2f AK %5.2f MP %.2f AS %6.2f/s IS %6.2f/s LI %f, DONE" % (chainid, ntest, nfail, nkept, nstay, IK, AK, MP, AS, IS, LI))
    return models, datas, weights, llks
#########################################################
if __name__ == "__main__":
    from srfpython.graphictools.gutils import *
    from srfpython.processingtools.multipro7 import *
    from utils import Timer
    from srfpython.depthinversion.utils import histogram2d

    #########################################################
    def test1():
        # test 1D pdfs
        v = np.linspace(-1., 1., 10000)

        if 0:
            # nans expected
            l0 = LogUni(0., 0.1, k=20., nanbehavior=2)
            l1 = LogGauss(0.05, 100000., 0., 0.1, k=20., nanbehavior=2)
            l2 = LogGauss(0., 0.1, 0., 0.1, k=20., nanbehavior=2)

            V = np.copy(v)
            V[-1] = np.nan
        else:
            # no nans expected
            l0 = LogUni(0., 0.1, k=20., nanbehavior=0)
            l1 = LogGauss(0.05, 100000., 0., 0.1, k=20., nanbehavior=0)
            l2 = LogGauss(0., 0.1, 0., 0.1, k=20., nanbehavior=0)
            V = np.copy(v)

        gcf().add_subplot(211)
        gca().plot(v, l0.calln(V), label="l0")
        gca().plot(v, l1.calln(V), label="l1")
        gca().plot(v, l2.calln(V), label="l2")
        gca().legend()
        gcf().add_subplot(212, sharex=gca())
        gca().plot(v, np.exp(l0.calln(V)), label="l0")
        gca().plot(v, np.exp(l1.calln(V)), label="l1")
        gca().plot(v, np.exp(l2.calln(V)), label="l2")
        gca().legend()
        showme(False)
    # ________________________________________________________
    def test2():
        # test ND pdfs
        x = np.linspace(0., 1., 200)
        y = np.linspace(0., 2., 215)
        X, Y = np.meshgrid(x, y)
        # l = LogUniND([0.1, 0.1], [0.3, 0.5], k = 20.)
        # l = LogGaussND([0.2, 0.2], [0.1, 1000.], [0.1, 0.1], [0.3, 0.5], k = 20.)
        l = LogGaussNDCov([0.2, 0.2], [0.1, 0.2], [0.1, 0.1], [0.3, 0.5], rho=[[1.0, 0.9], [0.9, 1.0]], k=20.,
                          nanbehavior=0)
        Z = l.callargs(X, Y)
        gcf().add_subplot(121)
        gca().pcolormesh(X, Y, Z)
        gcf().add_subplot(122, sharex=gca(), sharey=gca())
        gca().pcolormesh(X, Y, np.exp(Z))
        showme(False)

        #########################################################
    # ________________________________________________________
    def test211():
        # nans expected
        l = LogGaussND(
            vmeans=[0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                    0.2, 0.2, 0.2],
            vuncs=[0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1,
                   0.2, 0.1, 0.2],
            vinfs=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                   0.1, 0.1, 0.1],
            vsups=[0.3, 0.5, 0.3, 0.5, 0.3, 0.5, 0.3, 0.5, 0.3, 0.5, 0.3, 0.5, 0.3, 0.5, 0.3, 0.5, 0.3, 0.5, 0.3, 0.5, 0.3,
                   0.5, 0.3, 0.5],
            k=20.,
            nanbehavior=1)
        with Timer(''):
            for _ in xrange(1000):
                l.call1(
                    [np.nan, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                         0.2, 0.2, 0.2, 0.2])
    # ________________________________________________________
    def test212():
        # no nans expected
        l = LogGaussND(
            vmeans=[0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                    0.2, 0.2, 0.2],
            vuncs=[0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1,
                   0.2, 0.1, 0.2],
            vinfs=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                   0.1, 0.1, 0.1],
            vsups=[0.3, 0.5, 0.3, 0.5, 0.3, 0.5, 0.3, 0.5, 0.3, 0.5, 0.3, 0.5, 0.3, 0.5, 0.3, 0.5, 0.3, 0.5, 0.3, 0.5, 0.3,
                   0.5, 0.3, 0.5],
            k=20.,
            nanbehavior=0)

        for _ in xrange(1000):
            l.call1(
                [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.2, 0.2, 0.2, 0.2, 0.2])
    # ________________________________________________________
    def test30():
        # test Metropolis 1D
        # l = LogGauss(vmean = 0., vunc = 1.0, vinf = -0.1, vsup = 3.0, k = 1000.0, nanbehavior=  0) #not adapted for metropolis
        l = LogGaussND(vmeans=[0.], vuncs=[1.0], vinfs=[-0.1], vsups=[3.0], k=1000.0, nanbehavior=0)
        x = np.linspace(-10., 10., 1000)
        y = np.exp([l([xx]) for xx in x])

        A = np.sum(y * (x[1] - x[0]))
        gca().plot(x, y / A)
        showme(False)

        def G(model): return np.array([0.])

        def logRHOD(data): return 0.

        logRHOM = l
        M0 = np.array([0.])
        MSTD = np.array([3.])
        models, _, weights, _ = metropolis(M0, MSTD, G, 1, logRHOD, logRHOM, nkeep=100000, HL=1000, IK0=0.25,
                                           MPMIN=1.0, MPMAX=1.0)
        models = models.repeat(weights, axis=0)
        #p = PDF(models[:, 0], 100)
        hist, edges = np.histogram(models[:, 0], bins = 100, density = True)
        gca().plot(0.5 * (edges[:-1] + edges[1:]), hist)
        showme(False)
    # ________________________________________________________
    def test31():
        # test Metropolis 2D, parallel

        x = np.linspace(0.05, 0.3, 200)
        y = np.linspace(0.05, 0.4, 215)
        X, Y = np.meshgrid(x, y)
        # l = LogUniND([0.1, 0.1], [0.3, 0.5], k = 20.)
        # l = LogGaussND([0.2, 0.2], [0.025, 1000.], [0.1, 0.1], [0.3, 0.5], k = 20.)
        # l = LogGaussNDCov([0.2, 0.2], [0.01, 0.2], [0.1, 0.1], [0.2, 0.5], rho = [[1.0, 0.9], [0.9, 1.0]], k = 20.)
        l = LogGaussNDCov([0.2, 0.2], [0.1, 0.2], [0.1, 0.1], [0.2, 0.5], rho=[[1.0, 0.9], [0.9, 1.0]], k=20.,
                          nanbehavior=0)
        Z = l.callargs(X, Y)
        ax1 = gcf().add_subplot(131)
        gca().pcolormesh(X, Y, Z, cmap=plt.cm.jet)#cmap2isocolor(Z, plt.cm.jet))
        ax2 = gcf().add_subplot(132, sharex=gca(), sharey=gca())
        gca().pcolormesh(X, Y, np.exp(Z))
        ax3 = gcf().add_subplot(133, sharex=gca(), sharey=gca());
        ax3 = gca()
        showme(False)

        def gen():
            for nchain in xrange(48):
                M0 = np.random.randn(2)
                MSTD = np.asarray([0.05, 0.05])
                yield Job(nchain, M0, MSTD, nkeep=1500)

        def fun(worker, chainid, M0, MSTD, nkeep):
            def G(model): return np.array([0.])

            def logRHOD(data): return 0.

            logRHOM = l
            models, _, weights, _ = metropolis(M0, MSTD, G, 1, logRHOD, logRHOM, nkeep=nkeep,
                                               normallaw=worker.randn, unilaw=worker.rand, chainid=chainid)
            models = models.repeat(weights, axis=0)

            return models[:, 0], models[:, 1]

        X, Y = [], []
        with MapAsync(fun, gen()) as ma:
            for _, (XX, YY), _, _ in ma:
                X = np.concatenate((X, XX))
                Y = np.concatenate((Y, YY))

        X, Y, H = histogram2d(X, Y, x[::2], y[::2])
        ax3.pcolormesh(X, Y, H)

        showme(False)
    # ________________________________________________________
    def test32():

        # gaussian regression example : assume incompatible data, and multiple solutions
        # -----------------------
        class Theo(object):
            def __init__(self, xpoints):
                self.xpoints = xpoints

            def __call__(self, m):
                # theory
                # d = m[0] * self.xpoints + m[1]
                d = np.exp(-0.5 * ((self.xpoints - m[0]) / m[1]) ** 2.)
                return d

        # -----------------------
        # data
        amin, amax = 0., 1.0
        bmin, bmax = 0.05, 0.2
        x = np.linspace(0., 1.0, 24)
        G = Theo(x)
        # -----------------------
        y = np.zeros_like(x)
        y[:12] = G([0.25, 0.05])[:12] + 0.05 * np.random.randn(len(x[:12]))
        y[12:] = G([0.75, 0.08])[12:] + 0.05 * np.random.randn(len(x[12:]))
        s = 0.3 * np.ones_like(x)
        plt.figure(figsize=(12, 6))
        ax1 = gcf().add_subplot(121, xlabel="mean", ylabel="std");
        ax2 = gcf().add_subplot(122, xlabel="x", ylabel="y");

        for xx, yy, ss in zip(x, y, s):
            ax2.plot([xx, xx], [yy - ss, yy + ss], "_-", color=[0.5, .5, .5], alpha=1.0)
        ax2.plot(x, y, "wo", markeredgecolor="k")
        showme(0)
        # -----------------------
        logRHOD = LogGaussND( \
            vmeans=y,
            vuncs=s,
            vinfs=y * 0. - 10.,
            vsups=y * 0. + 10.)
        logRHOM = LogUniND( \
            vinfs=[amin, bmin],
            vsups=[amax, bmax])
        # -----------------------
        if False:
            M0 = np.array([np.random.rand() * (amax - amin) + amin, np.random.rand() * (bmax - bmin) + bmin])
            MSTD = np.array([0.1, 0.1])
            nkeep = 10000
            models, datas, weights, llks = metropolis(M0, MSTD, G, len(G.xpoints), logRHOD, logRHOM, nkeep=nkeep)
            models = models.repeat(weights, axis=0)
            datas = datas.repeat(weights, axis=0)
            llks = llks.repeat(weights)
            A = models[:, 0]
            B = models[:, 1]
        else:  # parallel
            def gen():
                for nchain in xrange(12):
                    M0 = np.array(
                        [np.random.rand() * (amax - amin) + amin, np.random.rand() * (bmax - bmin) + bmin])
                    MSTD = np.array([0.1, 0.1])
                    nkeep = 10000
                    yield Job(nchain, M0, MSTD, nkeep=nkeep)

            def fun(worker, chainid, M0, MSTD, nkeep):
                G = Theo(x)
                logRHOD = LogGaussND( \
                    vmeans=y,
                    vuncs=s,
                    vinfs=y * 0. - 10.,
                    vsups=y * 0. + 10.)
                logRHOM = LogUniND( \
                    vinfs=[amin, bmin],
                    vsups=[amax, bmax])
                models, datas, weights, llks = metropolis(M0, MSTD, G, len(G.xpoints), logRHOD, logRHOM,
                                                          nkeep=nkeep, normallaw=worker.randn, unilaw=worker.rand,
                                                          chainid=chainid)
                models = models.repeat(weights, axis=0)
                datas = datas.repeat(weights, axis=0)
                llks = llks.repeat(weights)
                return models, datas, weights, llks

            models = None
            datas = None
            with MapAsync(fun, gen()) as ma:
                for _, (ms, ds, ws, ls), _, _ in ma:
                    if models is None:
                        models, datas = ms, ds
                    else:
                        models = np.concatenate((models, ms), axis=0)
                        datas = np.concatenate((datas, ds), axis=0)
            A = models[:, 0]
            B = models[:, 1]

        # ----------------------- model histogram
        X, Y, H = histogram2d(A, B, np.linspace(amin, amax, 100), np.linspace(bmin, bmax, 100))#, normalization="pdf")
        cmap = cmap=plt.cm.CMRmap #cmap2isocolor(H, cmap=plt.cm.CMRmap)
        ax1.pcolormesh(X, Y, H, cmap=cmap)
        showme(False)
        # ----------------------- data histogram
        xbins = np.linspace(x.min(), x.max(), 100)
        ybins = np.linspace((y - s).min(), (y + s).max(), 110)
        Xbins, Ybins = np.meshgrid(xbins, ybins)
        Z = np.zeros_like(Xbins)
        J = np.arange(len(xbins))
        G1 = Theo(xbins)
        for a, b in zip(A, B):
            y = G1([a, b])
            I = np.argmin(abs(y - Ybins), axis=0)
            Z[I, J] += 1
        cmap = plt.cm.CMRmap #cmap2isocolor(Z, cmap=plt.cm.CMRmap)
        ax2.pcolormesh(xbins, ybins, Z, cmap=cmap)

        showme(False)
    # ________________________________________________________
    def test33():
        # test for failing theory or nans in data array
        # use the Programming Error exception to test the behavior of the code
        # use failure test or happy failure to mimic unsuited models
        nofail = False  # if True, the chain will continue running also if programming error is suspected
        nanbehavior = 2  # recommended : means that models or datas with nans are penalized according to the number of nans

        # 0 means that nans are producing errors, so the metropolis algo cannot compare models since they are all equally bad
        # -----------------------
        class Theo(object):
            def __init__(self, xpoints):
                self.xpoints = xpoints

            def __call__(self, m):
                # theory
                if False:
                    raise Exception('Programming Error test')

                if m[1] >= 0.16:
                    raise Exception('Theory failure test')
                    return np.nan * np.zeros_like(self.xpoints)
                elif m[1] <= 0.08:
                    return np.nan * np.zeros_like(self.xpoints)
                elif np.random.rand() < 0.01:
                    raise Exception('Happy failure test : the failure that occurs just for fun')
                    return np.nan * np.zeros_like(self.xpoints)
                d = np.exp(-0.5 * ((self.xpoints - m[0]) / m[1]) ** 2.)
                return d

        # -----------------------
        # data
        amin, amax = 0., 1.0
        bmin, bmax = 0.05, 0.2
        x = np.linspace(0., 1.0, 12)
        G = Theo(x)
        # -----------------------
        y = np.zeros_like(x)
        y = np.exp(-0.5 * ((x - 0.5) / 0.1) ** 2.) + 0.1 * np.random.randn(len(x))
        s = 0.3 * np.ones_like(x)
        plt.figure(figsize=(12, 6))
        ax1 = gcf().add_subplot(121, xlabel="mean", ylabel="std");
        ax1 = gca()
        ax2 = gcf().add_subplot(122, xlabel="x", ylabel="y");
        ax2 = gca()

        for xx, yy, ss in zip(x, y, s):
            ax2.plot([xx, xx], [yy - ss, yy + ss], "_-", color=[0.5, .5, .5], alpha=1.0)
        ax2.plot(x, y, "wo", markeredgecolor="k")
        showme(0)
        # -----------------------
        if True:
            def gen():
                for nchain in xrange(12):
                    M0 = np.array(
                        [np.random.rand() * (amax - amin) + amin, np.random.rand() * (bmax - bmin) + bmin])
                    MSTD = np.array([0.1, 0.1])
                    nkeep = 1000
                    yield Job(nchain, M0, MSTD, nkeep=nkeep)

            def fun(worker, chainid, M0, MSTD, nkeep):
                G = Theo(x)
                logRHOD = LogGaussND( \
                    vmeans=y,
                    vuncs=s,
                    vinfs=y * 0. - 10.,
                    vsups=y * 0. + 10.,
                    nanbehavior=nanbehavior)
                logRHOM = LogUniND( \
                    vinfs=[amin, bmin],
                    vsups=[amax, bmax],
                    nanbehavior=nanbehavior)
                models, datas, weights, llks = metropolis(M0, MSTD, G, len(G.xpoints), logRHOD, logRHOM,
                                                          nkeep=nkeep, normallaw=worker.randn, unilaw=worker.rand,
                                                          chainid=chainid, nofail=nofail)
                models = models.repeat(weights, axis=0)
                datas = datas.repeat(weights, axis=0)
                llks = llks.repeat(weights)
                return models, datas, weights, llks

            models = None
            datas = None
            with MapAsync(fun, gen()) as ma:
                for _, (ms, ds, ws, ls), _, _ in ma:
                    if models is None:
                        models, datas = ms, ds
                    else:
                        models = np.concatenate((models, ms), axis=0)
                        datas = np.concatenate((datas, ds), axis=0)
            A = models[:, 0]
            B = models[:, 1]
        # ----------------------- data histogram
        xbins = np.linspace(x.min(), x.max(), 100)
        ybins = np.linspace((y - s).min(), (y + s).max(), 110)
        Xbins, Ybins = np.meshgrid(xbins, ybins)
        Z = np.zeros_like(Xbins)
        J = np.arange(len(xbins))
        # G1 = Theo(xbins)
        for a, b in zip(A, B):
            try:
                y = np.exp(-0.5 * ((xbins - a) / b) ** 2.)  # G1([a, b])#G1([a, b])
            except KeyboardInterrupt:
                raise
            except:
                continue
            I = np.argmin(abs(y - Ybins), axis=0)
            Z[I, J] += 1
        cmap = plt.cm.CMRmap# cmap2isocolor(Z, cmap=plt.cm.CMRmap)
        ax2.pcolormesh(xbins, ybins, Z, cmap=cmap)
        # ----------------------- model histogram
        X, Y, H = histogram2d(A, B, np.linspace(amin - 0.33, amax + 0.33, 100),
                              np.linspace(bmin - 0.08, bmax + 0.08, 100))#, normalization="pdf")
        cmap = plt.cm.CMRmap#cmap2isocolor(H, cmap=plt.cm.CMRmap)
        ax1.pcolormesh(X, Y, H, cmap=cmap)
        showme(False)
    # ________________________________________________________
    def test4():
        """add constraint on secondary parameters"""
        # prior constraint on x    2*x
        l = LogUniND(vinfs=[-10., 0.0],
                     vsups=[10.0, 5.0],  # => x should be only between 0. and 2.5
                     k=1000.0,
                     nanbehavior=0)

        def logRHOM(model):
            x = model[0]
            y = 2. * x
            return l([x, y])

        x = np.linspace(-10., 10., 1000)
        y = np.exp([logRHOM([xx]) for xx in x])

        A = np.sum(y * (x[1] - x[0]))
        gca().plot(x, y / A)
        showme(False)

        def G(model): return np.array([0.])

        def logRHOD(data): return 0.

        M0 = np.array([0.])
        MSTD = np.array([8.])
        models, _, weights, _ = metropolis(M0, MSTD, G, 1, logRHOD, logRHOM, nkeep=100000, HL=1000, IK0=0.25,
                                           MPMIN=0.001, MPMAX=1000.0, adjustspeed=0.1)
        models = models.repeat(weights, axis=0)
        #p = PDF(models[:, 0], 200)
        hist, edges = np.histogram(models[:, 0], bins=200, density=True)
        #p.plot(gca(), alpha=0.5)
        gca().plot(0.5 * (edges[:-1] + edges[1:]), hist)
        showme(False)
    # ________________________________________________________
    def test41():
        """add constraint on secondary parameters, 2D
           model = x, y
           target = (x ** 2. + y ** 2.) ** 0.5 = 1.0 +- 0.1
           constraint = abs(x - y) follows a uniform law between two values

           > a model is [x, y]
           > the data is G(model) is [d]
           > the prior pdf on the model is bypassed so that additional constraints are added to some relations
             between the model parameters

        """

        # data pdf
        def G(model):
            x, y = model
            d = np.sqrt(x ** 2. + y ** 2.)
            return np.array([d])

        logRHOD = LogGaussND(vmeans=[1.0],
                             vuncs=[0.1],
                             vinfs=[0.5],
                             vsups=[1.5],
                             k=1000.0,
                             nanbehavior=0)

        # model pdf
        # prior constraint on x      y    abs(x-y)
        l = LogUniND(vinfs=[-10., -10.0, 0.5],
                     vsups=[10., 10.0, 1.5],
                     k=1000.0,
                     nanbehavior=0)

        def logRHOM(model):
            x, y = model
            return l([x, y, np.abs(x - y)])

        M0 = np.array([0., 0.])
        MSTD = np.array([1., 1.])
        models, _, weights, _ = metropolis(M0, MSTD, G, 1, logRHOD, logRHOM, nkeep=100000, HL=1000, IK0=0.25,
                                           MPMIN=0.001, MPMAX=1000.0, adjustspeed=0.1, debug=False)
        models = models.repeat(weights, axis=0)

        X, Y, H = histogram2d(models[:, 0], models[:, 1], 100, 100)#, normalization="pdf")
        cmap = plt.cm.CMRmap #cmap2isocolor(H, cmap=plt.cm.CMRmap)
        gca().pcolormesh(X, Y, H, cmap=cmap)
        gca().set_aspect(1.0)
        showme(False)
    #########################################################
    test41()
    pause()

