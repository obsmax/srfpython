import numpy as np
from numpy import log
import time


# ########################################################
class LogUni(object):
    """logarithm of a 1D, uniform, non normalized pdf (pseudo-uniform)
        I nead a pdf approaching a uniform distribution
        I want "bad" models to keep comparable, so I add a penalty that increases with increasing distance to the target
        zone

        f(x) =  | 0                          if     vinf <= x <= vsup
                | -1/2 * (x - vinf) / vstd   if     x < vinf
                | -1/2 * (x - vsup) / vstd   if     x > vsup
                | nanpenalty  or exception   if     x is nan or inf
                where vstd is (vsup - vinf) / k, k > 1
    """
    # ------------------------------------------------------
    def __init__(self, vinf, vsup, k=1000.0, nanbehavior=0):
        """
        :param vinf: float, lowest bound for uni pdf
        :param vsup: float, highest bound for uni pdf
        :param k: float, sharpness coefficients for the edges of the non-zero probability area
                  edges are defined as gaussian functions with std like (vsup - vinf) / k, k > 1
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
        self.k = float(k)
        self.vmean = 0.5 * (self.vinf + self.vsup)
        self.std = (self.vsup - self.vinf) / self.k
        self.nanbehavior = nanbehavior

        if self.nanbehavior == 0:
            self.raiseifnan = True

        elif self.nanbehavior == 1:
            # case 1: return a constant penalty, for ND pdfs, the penalty will be proportionnal to the number of nans
            self.raiseifnan = False
            self.nanpenalty = self.call1(self.vinf - 0.1 * (self.vsup - self.vinf)) #take the pdf value away from the 0 area

        elif self.nanbehavior == 2:
            # case 2: return severe penalty
            self.raiseifnan = False
            self.nanpenalty = self.call1(self.vinf - 10000. * (self.vsup - self.vinf))  # take the pdf value away from the 0 area
        else: raise ValueError('no such nan behavior implemented')

    # ------------------------------------------------------
    def __str__(self):
        return "{name}(vinf={vinf}, vsup={vsup}, k={k}, nanbehavior={nanbehavior})".format(
            name=self.__class__.__name__,
            vinf=self.vinf,
            vsup=self.vsup,
            k=self.k,
            nanbehavior=self.nanbehavior)

    # ------------------------------------------------------
    def calln(self, v):
        """v is an array with undefined length"""
        y = np.zeros_like(v)
        # -----------------------
        GOOD = ~(np.isnan(v) | np.isinf(v))
        if not GOOD.all():
            if self.raiseifnan: raise Exception('pdf got innapropriate values')
            else: y[~GOOD] = self.nanpenalty

        # -----------------------
        I = J = np.zeros_like(GOOD)

        I[GOOD] = (v[GOOD] < self.vinf)
        y[I] += -0.5 * ((v[I] - self.vinf) / self.std) ** 2.

        J[GOOD] = (self.vsup < v[GOOD])
        y[J] += -0.5 * ((v[J] - self.vsup) / self.std) ** 2.

        # -----------------------
        return y

    # ------------------------------------------------------
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

    # ------------------------------------------------------
    def __call__(self, v):
        """v is a scalar, see call1"""
        return self.call1(v)


# ########################################################
class LogGauss(object):
    """1D truncated gaussian, take the product of a gaussian and a pseudo-uniform distribution (see LogUni)
        the nan or inf behavior is handled by the pseudo uniform pdf

        f(x, *args) = LogUni(x, *args)  -1/2 * ((x - vmean) / vunc) ** 2.
    """
    def __init__(self, vmean, vunc, vinf, vsup, k = 1000., nanbehavior = 0):
        """
        :param vmean: float, mean value of the gaussian term, must be between vinf and vsup
        :param vunc: float, standard deviation of the gaussian term
        :param vinf: float, lower boundary of the uniform term
        :param vsup: float, upper boundary  of the uniform term
        :param k: float, see LogUni
        :param nanbehavior: int, see LogUni
        """
        assert not np.any(np.isnan([vmean, vunc, vinf, vsup, k]))
        assert not np.any(np.isinf([vmean, vunc, vinf, vsup, k]))
        assert vinf <= vmean <= vsup
        assert k > 0.
        self.vinf, self.vsup = vinf, vsup
        self.k, self.nanbehavior = float(k), nanbehavior
        self.vmean, self.vunc = vmean, vunc
        self.luni = LogUni(vinf, vsup, k, nanbehavior)

    # ------------------------------------------------------
    def __str__(self):
        return "{name}(vmean={vmean}, vunc={vunc}, vinf={vinf}, vsup={vsup}, k={k}, nanbehavior={nanbehavior})".format(
            name=self.__class__.__name__,
            vmean=self.vmean,
            vunc=self.vunc,
            vinf=self.vinf,
            vsup=self.vsup,
            k=self.k,
            nanbehavior=self.nanbehavior)

    # ------------------------------------------------------
    def calln(self, v):
        """v is an array with undefined length"""
        y = self.luni.calln(v) #handles nans or infs according to nanbehavior

        K = ~(np.isnan(v) | np.isinf(v))
        y[K] += -0.5 * ((v[K] - self.vmean) / self.vunc) ** 2.

        return y

    # ------------------------------------------------------
    def call1(self, v):
        """v is a scalar"""
        y = self.luni.call1(v)
        if np.isnan(v) or np.isinf(v): return y
        return y + -0.5 * ((v - self.vmean) / self.vunc) ** 2.

    # ------------------------------------------------------
    def __call__(self, v):
        return self.call1(v)


# ########################################################
class LogUniND(object):
    """N uniform laws, assume independent variables, sum of LogUni objects"""
    def __init__(self, vinfs, vsups, k = 1000., nanbehavior = 0):
        """
        :param vinfs: list of floats (one per dimension), see LogUni
        :param vsups: list of floats (one per dimension), see LogUni
        :param k: float, see LogUni
        :param nanbehavior: int, see LogUni
        """
        vinfs, vsups = [np.asarray(_, float) for _ in [vinfs, vsups]]
        assert len(vinfs) == len(vsups)
        self.N = len(vinfs)
        self.ls = [LogUni(vinf, vsup, k, nanbehavior) for vinf, vsup in zip(vinfs, vsups)]

    # ------------------------------------------------------
    def __str__(self):
        return "{name}(vinfs={vinfs}, vsups={vsups}, k={k}, nanbehavior={nanbehavior})".format(
            name=self.__class__.__name__,
            vmean=self.vmean,
            vunc=self.vunc,
            vinf=self.vinf,
            vsup=self.vsup,
            k=self.k,
            nanbehavior=self.nanbehavior)

    # ------------------------------------------------------
    def callpoints(self, points):
        """points = 2D array : raws = list of points, columns = dimensions"""
        y = np.zeros_like(points.shape[0], float)
        for n, l in enumerate(self.ls):
            y += l(points[:, n])
        return y

    # ------------------------------------------------------
    def callargs(self, *args):
        """args is like (X0, X1, X2, ..., XM-1) where Xi are arrays with same shapes"""
        assert len(args) == self.N
        for arg in args[1:]:
            assert arg.shape == args[0].shape

        y = np.zeros_like(arg)
        for n, (l, arg) in enumerate(zip(self.ls, args)):
            y += l.calln(arg)
        return y

    # ------------------------------------------------------
    def call1(self, model):
        #model like [x0, x1, x2, ..., xM-1] where xi are scalars
        assert len(model) == self.N
        y = sum([l(model[j]) for j, l in enumerate(self.ls)])
        return y

    # ------------------------------------------------------
    def __call__(self, model):
        #model like [x0, x1, x2, ..., xM-1] where xi are scalars
        return self.call1(model)


# ########################################################
class LogGaussND(LogUniND):
    """ND truncated gaussians, no covariance, sum of LogGauss objects"""
    def __init__(self, vmeans, vuncs, vinfs, vsups, k = 1000., nanbehavior = 0):
        """
        :param vmeans: list of floats, see LogGauss
        :param vuncs: list of floats, see LogGauss
        :param vinfs: list of floats, see LogUni
        :param vsups: list of floats, see LogUni
        :param k: float, see LogUni
        :param nanbehavior: integer, see LogUni
        """
        vmeans, vuncs, vinfs, vsups = [np.asarray(_, float) for _ in [vmeans, vuncs, vinfs, vsups]]
        assert len(vinfs) == len(vsups) == len(vmeans) == len(vuncs)
        self.N = len(vinfs)
        self.ls = [LogGauss(vmean, vunc, vinf, vsup, k, nanbehavior) \
                        for vmean, vunc, vinf, vsup in zip(vmeans, vuncs, vinfs, vsups)]


# ########################################################
class LogGaussNDCov(LogUniND):
    """ND truncated gaussians, with covariance"""
    def __init__(self, vmeans, vuncs, vinfs, vsups, rho, k = 1000., nanbehavior = 0):
        """
        :param vmeans: list of floats, see LogGauss
        :param vuncs: list of floats, see LogGauss
        :param vinfs: list of floats, see LogUni
        :param vsups: list of floats, see LogUni
        :param rho: 2D array, linear correlation coefficient, between -1 and 1
        :param k: float, see LogUni
        :param nanbehavior: integer, see LogUni

        """
        vmeans, vuncs, vinfs, vsups, rho = [np.asarray(_, float) for _ in [vmeans, vuncs, vinfs, vsups, rho]]
        assert len(vinfs) == len(vsups) == len(vmeans) == len(vuncs)
        assert rho.shape == (len(vinfs), len(vinfs))
        assert np.all((-1.0 <= rho) & (rho <= 1.0))
        self.N = len(vinfs)

        self.luni = LogUniND(vinfs, vsups, k, nanbehavior)
        sX, sY = np.meshgrid(vuncs, vuncs)
        self.CDinv = np.linalg.inv(rho * sX * sY)
        self.vmeans = np.asarray(vmeans)

    # ------------------------------------------------------
    def callargs(self, *args):
        #args like (X0, X1, X2, ..., XM-1) where Xi are arrays with same shapes
        y = self.luni.callargs(*args) #includes penalty for nans or infs
        for n, model in enumerate(zip(*[arg.flat[:] for arg in args])):
            y.flat[n] += self.call1(model)
        return y

    # ------------------------------------------------------
    def call1(self, model):
        #model like [x0, x1, x2, ..., xM-1] where xi are scalars
        M = np.asarray(model) - self.vmeans
        y = -.5 * np.dot(np.dot(M, self.CDinv), M.T) + self.luni.call1(model)
        return y

    # ------------------------------------------------------
    def __call__(self, model):
        #model like [x0, x1, x2, ..., xM-1] where xi are scalars
        return self.call1(model)


# ########################################################
def metropolis(M0, MSTD, G, ND, logRHOD, logRHOM, nkeep,
               normallaw = np.random.randn, unilaw = np.random.rand,
               chainid=1, HL = 100, IK0 = 0.25,
               MPMIN = 1.e-6, MPMAX = 1e6, adjustspeed=0.01, nofail=False, debug=False, verbose=True):
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
    notations :
        IK = instantaneous keep ratio
        AK = average keep ratio
        IS = instantaneous speed
        AS = average speed
        MP = master proposal
        LI = log likelyhood
    #-------------------
    variables:
        ntest    = number of tested models since chain start
        nfail    = number of models for which data could not be computed since chain start (the chain still samples the prior pdf)
        nkept    = number of distinct(!) models retained since chain start
        nstay    = number of subsequent rejection (= current weight of the last generated model)
        icurrent = the position in array models of the last generated model

    """
    summary = """chain%5d %5s kept%5d/%5d fail%5d AK%5.2f MP%5.2f AS%7.2f/s LI %f"""
    # ----
    nfail = 0

    # ----
    def call(M, nfail):
        try:
            D = G(M)
            L = logRHOM(M) + logRHOD(D)
        except KeyboardInterrupt: raise
        except Exception:
            if debug:
                # helps user to understand why the theory fails
                raise

            #could not compute data array, still sample the prior pdf to converge toward successfull regions, add constant penalty
            nfail += 1
            D = np.zeros(ND, float) * np.nan
            L = logRHOM(M) + -1e20 #bad models remain comparable, try to guide the chain away from these zones based on prior info
        return D, L, nfail

    # ----
    MI = M0
    DI, LI, nfail = call(MI, nfail)
    # ----
    ntest = nkept = nstay = 0
    MP    = 1.0 #master proposal, used to scale the markov chain proposal stds
    Ikept = np.zeros(HL, bool) #chain status history
    start = time.time()
    lasttime, lastntest = start, 0

    # ---- chain data
    # store every distinct model and data of the chain
    models  = np.empty((nkeep+1, len(MSTD)), float)
    datas   = np.empty((nkeep+1, len(DI)), float)
    llks    = np.zeros(nkeep+1, float) + -np.inf
    weights = np.zeros(nkeep+1, int)

    icurrent            = 0
    models[icurrent, :] = MI
    datas[icurrent, :]  = DI
    llks[icurrent]      = LI
    weights[icurrent]  += 1
    # ----
    while nkept < nkeep:
        ntest += 1
        # ----------------------
        if not ntest % HL and ntest:
            # reevaluate stats and master proposal
            IK  = Ikept.sum() / float(HL)  # instantaneous keep ratio, gives the recent acceptance of the chain
            AK  = nkept / float(ntest)     # average keep ratio
            MP *= 1. - adjustspeed * (IK0 - IK)   # master proposal, used to scale the markov chain proposal stds
            MP = np.clip(MP, MPMIN, MPMAX)
            t = time.time()
            AS  = ntest / (t - start)                  # Average test speed
            IS  = (ntest - lastntest) / (t - lasttime) # Instantaneous test speed
            lasttime, lastntest = t, ntest

        # ----------------------
        if not ntest % (10 * HL) and ntest:
            #reevaluate proposal stds
            I = np.arange(np.max([0, icurrent - 10 * HL]), icurrent + 1)
            MSTD = np.std(models[I, :].repeat(weights[I], axis = 0), axis = 0)

        # ----------------------
        if nfail >= nkeep and nfail >= ntest and not nofail:
            #all attempts to call G failed, it might be due to a programming error...
            print (summary % (chainid, "ERROR", nkept, ntest, nfail, AK, MP, AS, LI)) #PRESUMED ERROR IN THEORY
            while icurrent:
                G(models[icurrent, :]) #run G to reproduce the error message
                icurrent -= 1
            raise Exception('should never occur')

        # ----------------------
        if nstay >= nkeep:
            #stuck signal
            print (summary % (chainid, "STUCK", nkept, ntest, nfail, AK, MP, AS, LI))
            return models[:icurrent, :], datas[:icurrent, :], weights[:icurrent], llks[:icurrent]

        # ----------------------
        if not ntest % (10 * HL) and ntest: #True
            #report stats
            if verbose:
                print \
                ("chainid %5d ntest %5d nfail %5d nkept %5d nstay %5d IK %5.2f AK %5.2f MP %.2f AS %6.2f/s IS %6.2f/s LI %f" % \
                (chainid, ntest, nfail, nkept, nstay, IK, AK, MP, AS, IS, LI))
            #print (MP, " ".join(['%.4f' % _ for _ in MSTD / MSTD.max()]))

        # ----------------------
        MII = MI + MP * MSTD * normallaw(len(M0))
        DII, LII, nfail = call(MII, nfail)

        # ----------------------
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

        # ----------------------
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

    print (summary % (chainid, "DONE", nkept, ntest, nfail, AK, MP, AS, LI))
    return models, datas, weights, llks


#
