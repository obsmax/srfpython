#!/usr/bin/env python

from matplotlib import pyplot as plt
from matplotlib import colors
import numpy as np


# -------------------------------
# UTILS
# -------------------------------
def array2cmap(X):
    """convert an array X with RGB data into a colormap"""
    N = X.shape[0]
    
    r = np.linspace(0., 1., N+1)
    r = np.sort(np.concatenate((r, r)))[1:-1]

    rd = np.concatenate([[X[i, 0], X[i, 0]] for i in xrange(N)])
    gr = np.concatenate([[X[i, 1], X[i, 1]] for i in xrange(N)])
    bl = np.concatenate([[X[i, 2], X[i, 2]] for i in xrange(N)])
    
    rd = tuple([(r[i], rd[i], rd[i]) for i in xrange(2 * N)])
    gr = tuple([(r[i], gr[i], gr[i]) for i in xrange(2 * N)])
    bl = tuple([(r[i], bl[i], bl[i]) for i in xrange(2 * N)])

   
    cdict = {'red': rd, 'green': gr, 'blue': bl}
    return colors.LinearSegmentedColormap('my_colormap', cdict, N)


# -------------------------------
def cmapA2B(A = [0.5, 0.5, 0.5], B = [1.0, 0., 0.], N = 256):
    """create linear colorbar from colors A(r,g,b) to B(r,g,b)"""
    X = np.zeros((N, 3), float)
    for i in xrange(3):
        X[:, i] = np.linspace(A[i], B[i], N)
    return array2cmap(X)


# -------------------------------
def display_cmap(ax, cmap):
    X = cmap(np.arange(256))
    r, g, b = X[:, 0], X[:, 1], X[:, 2]
    ax.plot(r, color="r")
    ax.plot(g, color="g")
    ax.plot(b, color="b")
    ax.set_xlim(0, 255)
    cb = plt.cm.ScalarMappable(norm=None, cmap=cmap)
    cb.set_array([0, 255])
    cb.set_clim((0, 255))
    ax.figure.colorbar(cb, orientation="horizontal")


# -------------------------------
def stretchcmap(cmap, x, y=None, show=True):
    """modifies the repartition of colors along the colorbar using checkpoints
    input :
        cmap = a colormab object, e.g. plt.get_cmap('jet')
        x    = iterable, sorted, 0 to 1, horizontal positions of checkpoints
        y    = None or iterable, same len as x, 0 to 1, vertical positions of checkpoints
        show = bool, activate display for qc
    output :
        cmapout = new colormap

    >> demo : try
        cmap1 = stretchcmap(plt.get_cmap("spectral"), [0., 0.1, 1.0], [0., 0.4, 1.0], show = True)
    """

    x = np.asarray(x, float)
    assert x[0] == 0. and x[-1] == 1. and len(x) >= 3
    assert np.all(x[1:] >= x[:-1])
    if y is None:
        # default behavior
        y = np.linspace(0., 1., len(x))
    else:
        y = np.asarray(y, float)
        assert len(y) == len(x)
        assert y[0] == 0.
        assert y[-1] == 1.

    # ------------------------
    i = np.arange(256)
    X = cmap(i)
    r, g, b = X[:, 0], X[:, 1], X[:, 2]
    j = np.interp(xp=x * 255, fp=y * 255, x=i)
    r = np.interp(xp=i, fp=r, x=j)
    g = np.interp(xp=i, fp=g, x=j)
    b = np.interp(xp=i, fp=b, x=j)
    X[:, 0], X[:, 1], X[:, 2] = r, g, b
    cmapout = array2cmap(X)
    # ------------------------
    if show:
        fig = plt.figure(figsize=(6, 6))
        fig.subplots_adjust(left=0.2, bottom=0.2)
        ax = fig.add_subplot(111)
        pos = ax.get_position()
        caxv = ax.figure.add_axes((pos.x0 - 0.1 - 0.05 * pos.width, pos.y0, 0.05, pos.height))
        caxh = ax.figure.add_axes((pos.x0, pos.y0 - 0.1 - 0.05 * pos.height, pos.width, 0.05))
        ax.plot(x, y, "ko-")
        ax.plot([0., 1.], [0., 1.], "k--")

        cbv = plt.cm.ScalarMappable(norm=None, cmap=cmap)
        cbv.set_array([0, 255])
        cbv.set_clim((0, 255))
        fig.colorbar(cbv, cax=caxv, ticks=[])

        cbh = plt.cm.ScalarMappable(norm=None, cmap=cmapout)
        cbh.set_array([0, 255])
        cbh.set_clim((0, 255))
        fig.colorbar(cbh, cax=caxh, ticks=[], orientation="horizontal")
        ax.grid(True)
    # ------------------------
    return cmapout


# -------------------------------
# CMAPS
# -------------------------------
# def randcmap(N=256):
#     def randcolor():
#         r, g, b = np.random.rand(), np.random.rand(), np.random.rand()
#         return r, g, b
#
#     def getcolors(N):
#         r, g, b = randcolor()
#         colors = [randcolor()]
#         while len(colors) < N:
#             r, g, b = randcolor()
#             keepcolor = True
#             for i in xrange(len(colors)):
#                 ri, gi, bi = colors[i]
#                 if ((r - ri) ** 2. + (g - gi) ** 2. + (b - bi) ** 2.) ** .5 < .2:
#                     keepcolor=False
#                     break
#             if keepcolor:
#                 colors.append((r, g, b))
#         return colors
#
#     X = []
#     for r,g,b in getcolors(N):
#         X.append([r, g, b])
#
#     return array2cmap(np.array(X))


# -------------------------------
def linecmap(N=12):
    def linecolors():
        mycolors = np.array([\
            [1.,    0.,     0.],#r
            [0.,    1.,     0.],#g      
            [0.,    0.,     1.],#b
            [1.,    1.,     0.],#y       
            [0.,    1.,     1.],#c
            [1.,    0.,     1.],#p
            [1.,    .5,     0.],#orange
            [.8,    0.,     .5],#purple
            [.7,    1.,     0.],
            [0.,    1.,     .5],
            [.5,    0.,     1.],
            [0.,    .5,     1.],
            [1.,    .5,     .5],
            [.5,    1.,     .5],
            [.5,    .5,     1.] ])

        for n in xrange(mycolors.shape[0]):
            r, g, b = mycolors[n, :]
            yield r, g, b
        for n in xrange(mycolors.shape[0]):
            r, g, b = mycolors[n, :]
            yield r/1.5, g/1.5, b/1.5
        for n in xrange(mycolors.shape[0]):
            r, g, b = mycolors[n, :]
            yield r/3., g/3., b/3.
        yield 0., 0., 0.
        yield .3, .3, .3
        yield .6, .6, .6

    X = []
    gen = linecolors()
    for i in xrange(N):
        X.append(gen.next())

    X = np.array(X)
    return array2cmap(X)
###########################################################################################
def jetwk():
    cmap = plt.cm.get_cmap('jet')
    cmap.set_under('w')
    cmap.set_over('k')
    return cmap

import warnings
def yarg():
    warnings.warn('obsolet use plt.get_cmap')
    return plt.get_cmap('gray_r')


def toh():
    warnings.warn('obsolet use plt.get_cmap')
    return plt.get_cmap('hot_r')


def pamRMC():
    warnings.warn('obsolet use plt.get_cmap')
    return plt.get_cmap('CMRmap_r')


def cimsies():
    warnings.warn('obsolet use plt.get_cmap')
    return plt.get_cmap("seismic_r")


def tej():
    warnings.warn('obsolet use plt.get_cmap')
    return plt.get_cmap('jet_r')


def lartceps():
    warnings.warn('obsolet use plt.get_cmap')
    return plt.get_cmap('nipy_spectral_r')


def racn_tsig():
    warnings.warn('obsolet use plt.get_cmap')
    return plt.get_cmap('gist_ncar_r')


def spectralwide():
    l1 = plt.get_cmap('nipy_spectral')(np.arange(246))
    l2 = plt.get_cmap('Set1')(np.arange(5, 256))#[::-1])
    l = np.concatenate((l1, l2))    
    return array2cmap(l)


def ediwlartceps():
    return array2cmap(spectralwide()(np.arange(246 + 256 - 5)[::-1]))


def tomocmap(w = .2, g = 0.5):
    """

    :param w:.2 #largeur de la bande grise
    :param g: 1. #paleur du gris
    :return:
    """

    cdict = {'red':    ((0.0,     0.0, 1.0),
                        (.5 - w,  1.0, 1.0),                                     
                        (0.5,     g,   g),
                        (.5 + w,  0.0, 0.0),
                        (1.0,     0.0, 0.0)),
              'green': ((0.0,     0.0, 0.0),
                        (.5 - w,  1.0, 1.0),
                        (0.5,     g,   g),
                        (.5 + w,  1.0, 1.0),
                        (1.0,     0.0, 0.0)),
              'blue':  ((0.0,     0.0, 0.0),
                        (.5 - w,  0.0, 0.0),
                        (0.5,     g,   g),
                        (.5 + w,  0.5, 0.5),
                        (1.0,     1.0, 0.0))}
    return colors.LinearSegmentedColormap('my_colormap', cdict, 256)


def tomocmap1(w=0.05, W=.30, Wsat = 0.48, sat=.2, Ncolor = 256, reverse = True):
        """
        w     = witdh of the green area near 0
        W    = width of the yellow to cian area
        Wsat = width of the darken zone to indicate saturation
        sat  = darkness of the saturated zones 0 means black
        """
        assert 0. <= w <= W <= Wsat <= 0.5
        cdict = {'red':    ((0.0,       0.0, 0.0),
                            (.5 - Wsat, 0.0, 0.0),
                            (.5 - W,    0.0, 0.0),  
                            (.5 - w,    0.0, 0.0),
                            (.5 + w,    1.0, 1.0),               
                            (.5 + W,    1.0, 1.0),
                            (.5 + Wsat, sat, sat),
                            (1.0,       sat, sat)),
                            
                  'green': ((0.0,       0.0, 0.0),
                            (.5 - Wsat, 0.0, 0.0),
                            (.5 - W,    0.0, 0.0),
                            (.5 - w,    1.0, 1.0),
                            (.5 + w,    1.0, 1.0),
                            (.5 + W,    0.0, 0.0),
                            (.5 + Wsat,   0., 0.),
                            (1.0,       0.0, 0.0)),
                            
                  'blue':  ((0.0,       sat, sat),
                            (.5 - Wsat, sat, sat),
                            (.5 - W,    1.0, 1.0),
                            (.5 - w,    1.0, 1.0),
                            (.5 + w,    0.0, 0.0),               
                            (.5 + W,    0.0, 0.0),
                            (.5 + Wsat,   0., 0.),
                            (1.0,       0.0, 0.0))}
        cmap = colors.LinearSegmentedColormap('my_colormap', cdict, 256)
        if reverse: return array2cmap(cmap(np.arange(256)[::-1]))
        return cmap


def pamcomot():
    return array2cmap(tomocmap()(np.arange(256)[::-1]))


def cccfcmap(w=0.05, W=.4, g=.2, Ncolor = 256):
    cdict = {'red':    ((0.0,     0.0, 0.0),
                        (.5 - W,  0.0, 0.0),  
                        (.5 - w,    g,   g),
                        (.5 + w,  1-g, 1-g),               
                        (.5 + W,  1.0, 1.0),
                        (1.0,     1.0, 0.0)),
              'green': ((0.0,     0.0, 0.0),
                        (.5 - W,  0.0, 0.0),
                        (.5 - w,    g,   g),
                        (.5 + w,  1-g, 1-g),               
                        (.5 + W,  0.0, 0.0),
                        (1.0,     0.0, 0.0)),
              'blue':  ((0.0,     0.0, 1.0),
                        (.5 - W,  1.0, 1.0),
                        (.5 - w,    g,   g),  
                        (.5 + w,  1-g, 1-g),               
                        (.5 + W,  0.0, 0.0),
                        (1.0,     0.0, 0.0))}
    return colors.LinearSegmentedColormap('my_colormap', cdict, Ncolor)


def cccfcmap2(w=0.05, W=.25, g=.2, Ncolor = 256):
    cdict = {'red':    ((0.0,     0.0, 0.0),
                        (.5 - W,  0.0, 0.0),  
                        (.5 - w,    g,   g),
                        (.5 + w,  1-g, 1-g),               
                        (.5 + W,  1.0, 1.0),
                        (1.0,     1.0, 1.0)),
              'green': ((0.0,     0.0, 1.0),
                        (.5 - W,  0.0, 0.0),
                        (.5 - w,    g,   g),
                        (.5 + w,  1-g, 1-g),               
                        (.5 + W,  0.0, 0.0),
                        (1.0,     1.0, 1.0)),
              'blue':  ((0.0,     0.0, 1.0),
                        (.5 - W,  1.0, 1.0),
                        (.5 - w,    g,   g),  
                        (.5 + w,  1-g, 1-g),               
                        (.5 + W,  0.0, 0.0),
                        (1.0,     0.0, 0.0))}
    return colors.LinearSegmentedColormap('my_colormap', cdict, Ncolor)


def cccfcmap3(w=0.02, W=.25, g=.0, Ncolor = 256):
    cdict = {'red':    ((0.0,     0.0, 0.0),
                        (.5 - W,  0.0, 0.0),  
                        (.5 - w,    g,   g),
                        (.5 + w,  1-g, 1-g),               
                        (.5 + W,  1.0, 1.0),
                        (1.0,     1.0, 1.0)),
              'green': ((0.0,     0.0, 0.0),
                        (.5 - W,  0.0, 0.0),
                        (.5 - w,  1.0, 1.0),
                        (.5 + w,  1.0, 1.0),               
                        (.5 + W,  0.0, 0.0),
                        (1.0,     0.0, 0.0)),
              'blue':  ((0.0,     1.0, 1.0),
                        (.5 - W,  1.0, 1.0),
                        (.5 - w,  1-g, 1-g),  
                        (.5 + w,    g,   g),               
                        (.5 + W,  0.0, 0.0),
                        (1.0,     0.0, 0.0))}
    return colors.LinearSegmentedColormap('my_colormap', cdict, 256)


def cccfcmap4(w=0.02, W=.15, Wsat = 0.45, g=0., Ncolor = 256, reverse = False):
    cdict = {'red':    ((0.0,     0.0, 0.0),
                        (.5 - Wsat, 0., 0.),
                        (.5 - W,  0.0, 0.0),  
                        (.5 - w,    g,   g),
                        (.5 + w,  1-g, 1-g),               
                        (.5 + W,  1.0, 1.0),
                        (.5 + Wsat, 1., 1.),
                        (1.0,     1.0, 1.0)),
              'green': ((0.0,     0.0, 0.0),
                        (.5 - Wsat, 0., 0.),
                        (.5 - W,  0.0, 0.0),
                        (.5 - w,  1.0, 1.0),
                        (.5 + w,  1.0, 1.0),               
                        (.5 + W,  0.0, 0.0),
                        (.5 + Wsat, 1., 1.),
                        (1.0,     1.0, 1.0)),
              'blue':  ((0.0,     0.0, 0.0),
                        (.5 - Wsat, 0., 0.),
                        (.5 - W,  1.0, 1.0),
                        (.5 - w,  1-g, 1-g),  
                        (.5 + w,    g,   g),               
                        (.5 + W,  0.0, 0.0),
                        (.5 + Wsat, 1., 1.),
                        (1.0,     1.0, 1.0))}
    cmap = colors.LinearSegmentedColormap('my_colormap', cdict, 256)
    if reverse:
        return array2cmap(cmap(np.arange(256)[::-1]))
    return cmap


def rdcmap():
    cdict = {'red':    ((0.0,       0.5, 0.5),
                        (0.05,       1.0, 1.0),
                        (0.75,      1.0, 1.0),
                        (1.0,       1.0, 1.0)),
              'green': ((0.0,       0.5, 0.5),
                        (0.05,       0.0, 0.0),              
                        (0.75,      1.0, 1.0),
                        (1.0,       1.0, 1.0)),
              'blue':  ((0.0,       0.5, 0.5),    
                        (0.05,       0.0, 0.0),
                        (0.75,      0.0, 0.0),
                        (1.0,       1.0, 1.0))}

    return colors.LinearSegmentedColormap('my_colormap', cdict, 256)


def decorrcmap():
    cdict = {'red':    ((0.0,       0.0, 0.0),
                        (0.5,       1.0, 1.0),
                        (0.75,      1.0, 1.0),
                        (1.0,       1.0, 1.0)),
              'green': ((0.0,       0.0, 0.0),
                        (0.5,       0.0, 0.0),              
                        (0.75,      1.0, 1.0),
                        (1.0,       1.0, 1.0)),
              'blue':  ((0.0,       1.0, 1.0),    
                        (0.5,       0.0, 0.0),
                        (0.75,      0.0, 0.0),
                        (1.0,       1.0, 1.0))}
    return colors.LinearSegmentedColormap('my_colormap', cdict, 256)


def graysat(w=0.1):
    """saturated gray scale
    """
    d = ((0.0,       0.0,  0.0),
         (0.5 - w,   0.0,  0.0),
         (0.5 + w,   1.0,  1.0),
         (1.0,       1.0,  1.0))
    cdict = {'red': d, 'green': d, 'blue': d}
    return colors.LinearSegmentedColormap('my_colormap', cdict, 500)


def R_B(w = 0.7):
    cdict = {'red':    ((0.0,       0.0, 0.0),
                        (0.5,       w, 0.0),

                        (0.5,       0.0, 1.0),
                        (1.0,       1.0, 1.0)),

              'green': ((0.0,       0.0, 0.0),
                        (0.5,       w, 0.0),    
          
                        (0.5,       0.0, w),
                        (1.0,       0.0, 0.0)),

              'blue':  ((0.0,       1.0, 1.0),    
                        (0.5,       1.0, 0.0),

                        (0.5,       0.0, w),
                        (1.0,       0.0, 0.0))}
    return colors.LinearSegmentedColormap('my_colormap', cdict, 256)


def daylight(morning = 7., evening = 20.):
    m, e = morning / 24., evening / 24.


    cdict = {'red':    ((0.0,       0.0, 0.0),
                        (m,         1.0, 1.0),
                        (e,         1.0, 1.0),
                        (1.0,       0.0, 0.0)),

              'green': ((0.0,       0.0, 0.0),
                        (m,         1.0, 1.0),
                        (e,         1.0, 1.0),
                        (1.0,       0.0, 0.0)),

              'blue':  ((0.0,       0.0, 0.0),
                        (m,         1.0, 1.0),
                        (e,         1.0, 1.0),
                        (1.0,       0.0, 0.0))}
    return colors.LinearSegmentedColormap('my_colormap', cdict, 256)


def bazcmap():
    return None


def megawide():
   
    lsup = plt.get_cmap('nipy_spectral')(np.arange(10, 256))
    linf = plt.cm.jet(np.arange(230)[::4])
    linf = 0.33 * (linf - linf.min()) / (linf.max() - linf.min())#np.clip(linf - 0.4, 0., 1.)
    Ntra = 30
    ltra = np.meshgrid(np.zeros(4), np.linspace(0., 1., Ntra))[1]
    ltra = ltra * (lsup[0, :] - linf[-1, :])  +  linf[-1, :]
    l = np.concatenate((linf, ltra, lsup))

    ################"
    l = l[:-20, :]
    linf = l
    lsup = plt.cm.jet(np.arange(230)[::-4])
    lsup = 1. - 0.33 * (lsup - lsup.min()) / (lsup.max() - lsup.min())#np.clip(linf - 0.4, 0., 1.)
    Ntra = 30
    ltra = np.meshgrid(np.zeros(4), np.linspace(0., 1., Ntra))[1]
    ltra = ltra * (lsup[0, :] - linf[-1, :])  +  linf[-1, :]
    l = np.concatenate((linf, ltra, lsup))


    return array2cmap(l)


def gistncarb(reverse=False):
    X = plt.cm.gist_ncar(np.arange(256))
    X[:13,2] = 0.
    if reverse: X = X[::-1, :]
    return array2cmap(X)


def bracntsig():
    return gistncarb(reverse=True)


def test():
    X = plt.cm.gist_ncar(np.arange(256))
    X[:13,2] = 0.
    # X[200:,1], X[200:,2] = X[200:,2], X[200:,1]
    # if reverse: X = X[::-1, :]
    return array2cmap(X)


if __name__ == "__main__":

    # -------------------------------
    allcmaps = {"linecmap": linecmap(),
                "jetwk": jetwk(),
                "yarg": yarg(),
                "toh": toh(),
                "pamRMC": pamRMC(),
                "cimsies": cimsies(),
                "tej": tej(),
                "lartceps": lartceps(),
                "racn_tsig": racn_tsig(),
                "spectralwide": spectralwide(),
                "ediwlartceps": ediwlartceps(),
                "tomocmap": tomocmap(),
                "tomocmap1": tomocmap1(),
                "pamcomot": pamcomot(),
                "cccfcmap": cccfcmap(),
                "cccfcmap2": cccfcmap2(),
                "cccfcmap3": cccfcmap3(),
                "cccfcmap4": cccfcmap4(),
                "rdcmap": rdcmap(),
                "decorrcmap": decorrcmap(),
                "graysat": graysat(),
                "R_B": R_B(),
                "daylight": daylight(),
                "bazcmap": bazcmap(),
                "megawide": megawide(),
                "gistncarb": gistncarb(),
                "bracntsig": bracntsig(),
                "test": test()}

    from srfpython.standalone.display import *

    # ############ demo display_cmap
    plt.figure()
    display_cmap(gca(), gistncarb())
    
    # ############ demo stretchcmap
    cmap1 = stretchcmap(plt.get_cmap("spectral"), [0., 0.1, 0.9, 1.0], [0., 0.4, 0.75, 1.0], show = True)
    showme()

    # ############ dispaly all cmaps
    plt.figure()

    keys = allcmaps.keys()
    for n, (key, cmap) in enumerate(allcmaps.items()):
        gcf().add_subplot(1, len(keys), n + 1)
        cb = makecolorbar(cmap = cmap, vmin = 0., vmax = 1.)
        cax=gca()
        gcf().colorbar(cb, cax=cax, ticks = [])
        cax.set_title(key, horizontalalignment="left", verticalalignment="bottom", rotation=45)

    showme()







