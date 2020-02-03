import numpy as np

"""
function optimization using the Nelder-Mead simplex algorithm

inversion :
    + this algorithm can solve non linear inverse problems
      by maximizing the posterior pdf in a way similar
      to the steepest ascent method
    - it is highly sensitive to local maxima in the
      posterior pdf and is therefore poorly efficient for problems
      with multiple solutions if the user does not
      have a fairly good prior estimation of the solution
    + maximizing the posterior pdf allows us to use fancy
      prior pdfs (e.g. uniform distributions, truncated gaussians, ...)
      while linearized approaches (e.g. newton) which operate on
      a linearized version of the theory use gaussian pdfs
    + the convergence is relatively fast (few evaluations of the pdf)
    + the algorithm is less sensitive to small fluctuations of the pdf (noise)
      since the slope is estimated through simplexes that are more stable
      than numerical derivation
"""

# -----------------------------------------
def neldermead(M0, DM, G, ND, logRHOD, logRHOM,
               alpha=1.5,
               beta=0.5,
               gamma=2.0,
               maxiter=1000,
               interrupt=1e-3,
               verbose=True,
               debug=False):
    """
    :param M0: starting point, shape (NM, ) where NM is the number of dimensions in the model space
    :type M0:  numpy.ndarray
    :param DM: the length of the simplex edge when initializing
    :type DM: float
    :param G: theory function, callable : np.ndarray (model) -> np.ndarray (data)
    :type G: function
    :param ND: number of dimensions in the data space
    :type ND: integer
    :param logRHOD: data prior pdf, callable  : np.ndarray (data) -> float (log(pdf(data)))
    :type logRHOD: function
    :param logRHOM: model prior pdf, callable : np.ndarray (model) -> float (log(pdf(model)))
    :type logRHOM:  function
    :param alpha: reflection coefficient, a factor for moving the worst point of the simplex toward the maximum
    :type alpha: float
    :param beta: contraction coefficient, the coeff for moving backward
    :type beta:  float
    :param gamma: expansion  coefficient, if the slope is long enough, then the step is increased
    :type gamma: float
    :param maxiter: int, maximum number of iteration that can be done
    :type maxiter: int
    :param interrupt: float, if the relative improvement is below interrupt 10 times in a raw, then the inversion is interrupted
    :type interrupt: float
    :param verbose: to print message about the convergence
    :type verbose: bool
    :param debug: bool,  raise Exception in case of theory failure
    :type debug:
    :return models, datas, llks:
        models : 2D array with the best model in each simplex tested (1 per line)
                 all tested models are not returned
        datas  : 2D array with the data array corresponding to each line of models
        llks   : 1D array, the log likelihood values corresponding to each model in models
    """
    
    # ----
    def call(M):
        try:
            D = G(M)
            L = logRHOM(M) + logRHOD(D)
        except KeyboardInterrupt:
            raise
        except Exception:
            if debug:
                # helps user to understand why the theory fails
                raise

            # could not compute data array,
            # we return a constant penalty that is very high (-1e20)
            # to notice that this model is very bad
            # we still add a penalty for the prior model
            # so that bad models remain comparable
            D = np.zeros(ND, float) * np.nan # means : data could not be computed
            L = logRHOM(M) + -1e20 # constant penalty + prior penalty
        return D, -L # return -L to follow th original algo (minimization instead of maximization)

    # ----
    assert 0. < beta < 1.0 <= alpha < gamma

    Ndim = len(M0) # number of dimensions in the model space
    Npoints = Ndim + 1 # number of points in the simplex
    if hasattr(DM, "__iter__"):
        assert len(DM) == Ndim
        assert np.all(DM > 0.)
    else:
        DM *= np.ones(Ndim, float)

    # ----------------- first simplex
    # 1 point per column, 1 line per dimension (example in 3D)
    #      M0    DM     Mis
    # x   [a]   [dx]   [a   a+dx   a      a   ]
    # y   [b]   [dy]   [b   b      b+dy   b   ]
    # z   [c]   [dz]   [c   c      c      c+dz]
    Mis = np.zeros((Ndim, Npoints), float)
    Mis[:, 1:] = DM * np.eye(Ndim) # => put DM on the diagonal of a square array
    Mis = (Mis.T + M0).T

    # ----------------- evaluate first simplex
    Lis = np.zeros(Npoints, float) * np.nan # log likelihood at simplex points
    Dis = np.zeros((ND, Npoints), float) * np.nan # Data corresponding to Mis
    for j in xrange(Ndim + 1):
        # compute the data and likelihood of all the points in the simplex
        Dis[:, j], Lis[j] = call(Mis[:, j])

    # -----------------
    models  = np.zeros((maxiter+1, Ndim), float) * np.nan
    datas   = np.zeros((maxiter+1, ND), float) * np.nan
    llks    = np.zeros(maxiter+1, float) + -np.inf
    # -----------------
    nstay = 0
    best = -np.inf
    for niter in xrange(maxiter):

        h = np.argmax(Lis)  # worst point (index)
        l = np.argmin(Lis)  # best point (index)
        Pb = np.mean(Mis, axis=1)  # center of mass

        ###################
        # bestpoint = Mis[:, l]
        # bestvalue = -Lis[l] #- because we evaluated -llk
        # # yield the best point of the simplex
        # yield bestpoint, bestvalue, (Mis, Dis, Lis)  # use this to see the simplexes tested
        # simplices.append(np.concatenate((Mis.T, [Mis[:, 0]]), axis=0))
        models[niter, :] = Mis[:, l]
        datas[niter, :] = Dis[:, l]
        llks[niter] = -Lis[l]
        if verbose:
            print "neldermead : iter %6d llk %f " % (niter, llks[niter])
        ###################

        ###################
        if niter > 1 and abs((-Lis[l] - best) / best) <= interrupt:
            nstay += 1
        else:
            nstay = 0
        if nstay > 10:
            # no significant improvement since 10 iterations
            break
        best = -Lis[l]
        ###################

        # ## REFLECTION
        Mstar = (1. + alpha) * Pb - alpha * Mis[:, h]
        Dstar, Lstar = call(Mstar) #Lstar = -fun(Mstar)

        if Lstar < Lis[l]:
            # reflection has produced a new minimum (Mstar)
            # ## EXPANSION
            M2star = gamma * Mstar + (1. - gamma) * Pb
            D2star, L2star = call(M2star)  # L2star = -fun(M2star)  # M2star[0], M2star[1])

            if L2star < Lis[l]:
                # ## EXPANSION SUCCEDED
                Mis[:, h] = M2star
                Dis[:, h] = D2star
                Lis[h] = L2star
            else:
                # ## EXPANSION FAILED
                Mis[:, h] = Mstar
                Dis[:, h] = Dstar
                Lis[h] = Lstar
        else:
            if np.all(Lstar > np.concatenate((Lis[:h], Lis[h + 1:]))):
                # Lstar AND yh are still the two worst points of the simplex
                if Lstar < Lis[h]:
                    # yh was worse than Lstar, then take yh = min(yh, Lstar)
                    Mis[:, h] = Mstar
                    Dis[:, h] = Dstar
                    Lis[h] = Lstar
                    # ## CONTRACTION
                M2star = beta * Mstar + (1. - beta) * Pb
                D2star, L2star = call(M2star)  # M2star[0], M2star[1])

                if L2star > Lis[h]:
                    # ## CONTRACTION FAILED (the contracted point (L2star is worst than Lstar)
                    Mis = 0.5 * (Mis.T + Mis[:, l]).T
                else:
                    # ## CONTRACTION SUCCEDED
                    Mis[:, h] = M2star
                    Dis[:, h] = D2star
                    Lis[h] = L2star
            else:
                # reflection has simply produced a new point
                Mis[:, h] = Mstar
                Dis[:, h] = Dstar
                Lis[h] = Lstar

    models, datas, llks = models[:niter, :], datas[:niter, :], llks[:niter]
    # I = np.argsort(llks) #reorder by loglikelihood
    #return models[I, :], datas[I, :], llks[I]
    return models, datas, llks
