import numpy as np


# --------------------------------------
def gauss2D2(x, y, xo, yo, sigmaj, sigmin, alphadeg):
    """plus pratique, pas de coeff de correl
    xo, yo = max of the gaussian
    sigmaj = standard deviation in the major eigenvector direction
    sigmin = standard deviation in the minor eigenvector direction
    alphadeg = angle in degrees, counterclockwise from the positive x axis
    """
    assert sigmaj >= sigmin

    deg2rad = np.pi / 180.
    x = x - xo
    y = y - yo
    c =  np.cos(alphadeg * deg2rad)
    s =  np.sin(alphadeg * deg2rad)
    xp = (y * c - x * s) / sigmin
    yp = (x * c + y * s) / sigmaj
    return np.exp(-0.5 * (xp ** 2. + yp ** 2.))

# -----------------------------------------
class fun2D():
    def __init__(self, N=10):
        if True:
            self.N = N
            self.As = np.random.rand(self.N)
            self.xos, self.yos = zip(*[np.random.rand(2) * 20. - 10. for _ in xrange(self.N)])
            self.sigmamins, self.sigmajs = zip(*[np.sort(np.random.rand(2) * 1. + .5) for _ in xrange(self.N)])
            self.alphadegs = np.random.rand(self.N) * 360.
        else:
            self.N = 1
            self.As = np.array([1.])
            self.xos = np.array([+1.123456789])
            self.yos = np.array([-1.123456789])
            self.sigmamins = np.array([.1])
            self.sigmajs   = np.array([3.])
            self.alphadegs = np.array([45.])

    def __call__(self, m):
        x, y = m
        z = np.zeros_like(x)
        for n in xrange(self.N):
            z += self.As[n] * gauss2D2(x, y, self.xos[n], self.yos[n], self.sigmajs[n], self.sigmamins[n],
                                       self.alphadegs[n])

        return z


# -----------------------------------------
def NelderMead(fun, P0, dx,
               alpha=1.5, beta=0.5, gamma=2.0,
               maxiter=1000,
               interrupt=1e-3):
    """
    simplex method for maximization
    P0    = array, starting point
    dx    = float, the length of the simplex edge when initializing
    alpha = float, reflection coefficient, a factor for moving the worst point of the simplex toward the minimum
    gamma = float, expansion  coefficient, if the slope is long enough, then the step is increased
    beta  = float, contraction coefficient, the coeff for moving backward
    maxiter = int, maximum number of iteration that can be done
    interrupt = float, if the relative improvement is below interrupt 10 times in a raw, then the inversion is interrupted

    fun must be a callable such as Z = fun([X, Y]), Z a scallar
    here I call -fun instead of fun, because I want to maximize fun (original algorithm was for minimization, See Nelder and Mead 1965)
    """
    assert alpha >= 1.0
    assert gamma >= 1.0
    assert 0. < beta < 1.0

    Ndim = 2
    Npoints = 3
    # first simplex 1 point per column, 1 line per dimension
    Pis = np.zeros((Ndim, Npoints), float)
    for _ in xrange(Npoints):
        Pis[:, _] = np.array(P0)

    Pis[0, 1] += dx
    Pis[1, 2] += dx

    yis = np.asarray([-fun(Pis[:, j]) for j in xrange(Pis.shape[1])])

    niter = 0
    nstay = 0
    best = -np.inf
    while niter < maxiter:
        niter += 1

        h = np.argmax(yis)  # worst point (index)
        l = np.argmin(yis)  # best point (index)
        Pb = np.mean(Pis, axis=1)  # center of mass

        ###################
        xtri = np.concatenate((Pis[0, :], [Pis[0, 0]]))
        ytri = np.concatenate((Pis[1, :], [Pis[1, 0]]))
        yield Pis[:, l], -yis[l], (xtri, ytri)  # yield the best point of the simplex
        ###################

        ###################
        if niter > 1 and abs((-yis[l] - best) / best) <= interrupt:
            nstay += 1
        else:
            nstay = 0
        if nstay > 10:
            break
        best = -yis[l]
        ###################


        # ## REFLECTION
        Pstar = (1. + alpha) * Pb - alpha * Pis[:, h]
        ystar = -fun(Pstar)

        if ystar < yis[l]:
            # reflection has produced a new minimum (Pstar)
            # ## EXPANSION
            P2star = gamma * Pstar + (1. - gamma) * Pb
            y2star = -fun(P2star)  # P2star[0], P2star[1])

            if y2star < yis[l]:
                # ## EXPANSION SUCCEDED
                Pis[:, h] = P2star
                yis[h] = y2star
            else:
                # ## EXPANSION FAILED
                Pis[:, h] = Pstar
                yis[h] = ystar
        else:
            if np.all(ystar > np.concatenate((yis[:h], yis[h + 1:]))):
                # ystar AND yh are still the two worst points of the simplex
                if ystar < yis[h]:
                    # yh was worse than ystar, then take yh = min(yh, ystar)
                    Pis[:, h] = Pstar
                    yis[h] = ystar
                    # ## CONTRACTION
                P2star = beta * Pstar + (1. - beta) * Pb
                y2star = -fun(P2star)  # P2star[0], P2star[1])

                if y2star > yis[h]:
                    # ## CONTRACTION FAILED (the contracted point (y2star is worst than ystar)
                    Pis = 0.5 * (Pis.T + Pis[:, l]).T
                else:
                    # ## CONTRACTION SUCCEDED
                    Pis[:, h] = P2star
                    yis[h] = y2star
            else:
                # reflection has simply produced a new point
                Pis[:, h] = Pstar
                yis[h] = ystar

if __name__ == "__main__":
    from tetedenoeud import *
    import time

    x = np.linspace(-10., 10., 200)
    y = np.linspace(-10., 10., 200)
    X, Y = np.meshgrid(x, y)
    f = fun2D()
    Z = f(np.asarray([X, Y]))

    plt.figure(figsize=(12, 6))
    ax0 = gcf().add_subplot(121);
    ax0 = gca()
    ax1 = gcf().add_subplot(122)

    ax0.pcolormesh(X, Y, Z, vmin=Z.min(), vmax=Z.max())
    #ax0.figure.show()
    ax0.set_aspect(1.0)

    for n, (Pi, yi, (xtri, ytri)) in enumerate(NelderMead(f, np.asarray([0., 0.]), dx=2., interrupt=1e-2)):
        ax0.plot(xtri, ytri, "k")
        ax0.plot(Pi[0], Pi[1], "k.")
        ax1.plot(n, yi, "ko")
        #gcf().canvas.draw()
        #time.sleep(0.1)

    showme()


