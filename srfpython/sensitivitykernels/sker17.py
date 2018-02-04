# from Herrmann.dispersion import *
import numpy as np
from srfpython.Herrmann.Herrmann import dispersion, dispersion_1, Timer, groupbywtm, igroupbywtm
from srfpython.utils import minmax

"""
srfker17, Maximilien Lehujeur, 01/11/2017
module to compute finite difference surface wave sensitivity kernels 
see documentation in function sker17
use __main__ for demo

see also Herrmann.dispersion.dispersion
"""

#_____________________________________
def sker17(ztop, vp, vs, rh, \
    waves, types, modes, freqs,
    delta = 0.05, norm = True, 
    h = 0.005, dcl = 0.005, dcr = 0.005):
    """sker17 : compute finite difference sensitivity kernels for surface waves dispersion curves 
    input: 
        -> depth model
        ztop, vp, vs, rh  : lists or arrays, see dispersion

        -> required dispersion points
        waves, types, modes, freqs : lists or arrays, see dispersion

        -> sensitivity kernel computation
        delta = float, > 0, relative perturbation to apply to each parameter and each layer
        norm  = bool, If false, I compute the simple derivative of the dispersion velocity
                                relative to a given parameter (e.g. dU(T)/dVs(z) where T is a given period and z is a given depth)
                      If true, I compute the ratio between the relative input perturbation (e.g. dVs(z)/Vs(z))
                               and the output perturbation (e.g. dU(T)/U(T))

        -> Herrmann's parameters, see CPS documentation
        h, dcl, dcr = passed to dispersion

    output:
        -> yields a tuple (w, t, m, F, DVADZ, DVADA, DVADB, DVADR) for each wave, type and mode
        w      = string, wave letter (L = Love or R = Rayleigh)
        t      = string, type letter (C = phase or U = group)
        m      = int, mode number (0= fundamental)
        F      = array, 1D, frequency array in Hz
        DVADZ  = array, 2D, [normed] sensitivity kernel relative to top depth of each layer (lines) and frequency (columns)
        DVADA  = array, 2D, [normed] sensitivity kernel relative to Pwave velocity of each layer (lines) and frequency (columns)
        DVADB  = array, 2D, [normed] sensitivity kernel relative to Swave velocity of each layer (lines) and frequency (columns)
        DVADR  = array, 2D, [normed] sensitivity kernel relative to density of each layer (lines) and frequency (columns)                
                 note that these arrays might contain nans
    see also :
        sker17_1
        dispersion
    """
    

    waves, types, modes, freqs = [np.asarray(_) for _ in waves, types, modes, freqs]
    nlayer = len(ztop)
    H = np.asarray(ztop)
    H[:-1], H[-1] = H[1:] - H[:-1], np.inf #layer thickness in km

    model0 = np.concatenate((ztop, vp, vs, rh))
    values0 = dispersion(ztop, vp, vs, rh, \
                       waves, types, modes, freqs, \
                       h = h, dcl = dcl, dcr = dcr)

    IZ  = np.arange(nlayer)
    IVP = np.arange(nlayer, 2*nlayer)    
    IVS = np.arange(2*nlayer, 3*nlayer)        
    IRH = np.arange(3*nlayer, 4*nlayer)        
    DVADP = np.zeros((4 * nlayer, len(waves)), float) * np.nan

    for i in xrange(1, 4 * len(ztop)):
        modeli = model0.copy()
        modeli[i] *= (1. + delta)
        ztopi, vpi, vsi, rhi = modeli[IZ], modeli[IVP], modeli[IVS], modeli[IRH]
        try:
            valuesi = dispersion(ztopi, vpi, vsi, rhi, \
                           waves, types, modes, freqs, \
                           h = h, dcl = dcl, dcr = dcr)        
        except CPiSDomainError as err: 
            print ("error during gradient computation %s" % str(err))
            continue
        except: raise
        if norm:
            DVADP[i, :] = ((valuesi - values0) / values0) / \
                          ((modeli[i] - model0[i]) / model0[i])
        else:
            DVADP[i, :] = (valuesi - values0) / (modeli[i] - model0[i]) 

    for w, t, m, F, Iwtm in groupbywtm(waves, types, modes, freqs, np.arange(len(waves))):
        DVADZ = DVADP[IZ, :][:, Iwtm]
        DVADA = DVADP[IVP, :][:, Iwtm]
        DVADB = DVADP[IVS, :][:, Iwtm]
        DVADR = DVADP[IRH, :][:, Iwtm]
        DVADZ, DVADA, DVADB, DVADR = [np.ma.masked_where(np.isnan(_), _) for _ in [DVADZ, DVADA, DVADB, DVADR]]
        
        yield w, t, m, F, DVADZ, DVADA, DVADB, DVADR

#_____________________________________
def sker17_1(ztop, vp, vs, rh, \
    Waves, Types, Modes, Freqs, **kwargs):
    """sker17_1 : same as sker17 with slightely more convenient input (no need to repeat wave, type and mode)

    Waves is like ['L', 'L', 'R']
    Types is like ['C', 'C', 'U']
    Modes is like [ 0,   1,   0 ]
    Freqs is like [fLC0, fLC1, fRU0] where f??? are frequency numpy arrays or lists
    
    see sker17 for detailed input and output arguments
    """
    waves, types, modes, freqs = igroupbywtm(Waves, Types, Modes, Freqs)
    for tup in sker17(ztop, vp, vs, rh, \
            waves, types, modes, freqs, **kwargs):
        yield tup

#_____________________________________
if __name__ == "__main__":
    """ DEMO """
    import matplotlib.pyplot as plt

    ###depth model
    ztop = [0.00, 0.25, 0.45, 0.65, 0.85, 1.05, 1.53, 1.80] #km, top layer depth
    vp   = [1.85, 2.36, 2.63, 3.15, 3.71, 4.54, 5.48, 5.80] #km/s
    vs   = [0.86, 1.10, 1.24, 1.47, 1.73, 2.13, 3.13, 3.31] #km/s
    rh   = [2.47, 2.47, 2.47, 2.47, 2.47, 2.58, 2.58, 2.63] #g/cm3

    ###dipsersion parameters
    def f(): return np.logspace(np.log10(0.1), np.log10(5.5), 35)
    Waves = ['R', 'R', 'R', 'R', 'L', 'L', 'L', 'L']
    Types = ['U', 'U', 'C', 'C', 'U', 'U', 'C', 'C']
    Modes = [ 0 ,  1,   0,   1,   0,   1,   0,   1 ]
    Freqs = [ f(), f(), f(), f(), f(), f(), f(), f()]

    ###compute dispersion curves
    with Timer('dispersion'):
        out = list(dispersion_1(ztop, vp, vs, rh, Waves, Types, Modes, Freqs))

    for w, t, m, fs, us in out:
        plt.gca().loglog(1. / fs, us, '+-', label = "%s%s%d" % (w, t, m))
    plt.gca().set_xlabel('period (s)')
    plt.gca().set_ylabel('velocity (km/s)')    
    plt.gca().grid(True, which = "major")
    plt.gca().grid(True, which = "minor")    
    plt.legend()
    plt.gcf().show()


    ###sensitivity kernels
    norm = True
    fig = plt.figure()  
    fig.subplots_adjust(wspace = 0.3, hspace = 0.3)
    fig.show()
    for w, t, m, F, DVADZ, DVADA, DVADB, DVADR in sker17_1(ztop, vp, vs, rh, \
        Waves, Types, Modes, Freqs, 
        norm = norm,   delta = 0.01, 
        h = 0.005, dcl = 0.005, dcr = 0.005):
        
        fig.clf()
        #------        
        ilayer = np.arange(DVADZ.shape[0]+1)-0.5
        iF = np.concatenate(([F[0] * 0.95], F * 1.05))   
        fig.suptitle('%s%s%d' % (w, t, m))
        
        #------
        vmin, vmax, cmap = -1., 1., plt.cm.RdBu
        ax1 = fig.add_subplot(221)        
        plt.pcolormesh(1./ iF, ilayer, DVADZ, vmin = vmin, vmax = vmax, cmap = cmap)
        plt.colorbar()
        ax2 = fig.add_subplot(222, sharex = ax1, sharey = ax1)        
        plt.pcolormesh(1./ iF, ilayer, DVADA, vmin = vmin, vmax = vmax, cmap = cmap)
        plt.colorbar()        
        ax3 = fig.add_subplot(223, sharex = ax1, sharey = ax1)        
        plt.pcolormesh(1./ iF, ilayer, DVADB, vmin = vmin, vmax = vmax, cmap = cmap)
        plt.colorbar()            
        ax4 = fig.add_subplot(224, sharex = ax1, sharey = ax1)
        plt.pcolormesh(1./ iF, ilayer, DVADR, vmin = vmin, vmax = vmax, cmap = cmap)
        plt.colorbar()            
        #------
        for ax, p in zip([ax1, ax2, ax3, ax4], ["Z", "Vp", "Vs", "rho"]):
            ax.set_xlim(minmax(1./ iF))
            ax.set_ylim(minmax(ilayer))        
            ax.set_xscale('log')
            ax.set_ylabel('layer number')
            ax.set_xlabel('period (s)')
            if norm:
                ax.set_title(r'$ \frac{d%s/%s}{d%s/%s} $' % (t, t, p, p))        
            else:
                ax.set_title('d%s/d%s' % (t, p))        
                
            if not ax.yaxis_inverted(): ax.invert_yaxis()            
        #------        
        fig.canvas.draw()
        raw_input('pause : press enter to plot the next wave type and mode')
    #--------------------    
    raw_input('bye')
    
