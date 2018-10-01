import matplotlib.pyplot as plt
from matplotlib import _pylab_helpers, ticker, axes, colors, __version__ as matplotlibversion
from matplotlib.offsetbox import AnchoredText
import numpy as np

# ################################################
gcf = plt.gcf
gca = plt.gca


# ______________________________________________
def pause(message="pause : press enter"):
    raw_input(message)


# ______________________________________________
def gcfs():
    "get current figures"
    return [w.canvas.figure for w in _pylab_helpers.Gcf.get_all_fig_managers()]


# ______________________________________________
def showme(block=True):
    figures = gcfs()
    for figure in figures:
        figure.canvas.draw()
        figure.show()
    if block: pause()


# ______________________________________________
def chftsz(obj, fontsize):
    "recursively change fontsize of an object including all its childrens"
    if hasattr(obj, "get_children"):
        for subobj in obj.get_children():
            chftsz(subobj, fontsize)
    if hasattr(obj, "set_fontsize"):
        obj.set_fontsize(fontsize)


# ################################################ tickers
def logtick(ax, axis='x', grid = True, color = "k", subs = [1., 2., 5.]):
    if matplotlibversion < "1.2":
        print "ERROR : logtick doesn't work with maptplotlib < 1.2"
        return

    def myfuncformatter(tick, tickposition=None):
        s = "%.6f" % tick
        if "." in s: s = s.rstrip('0').rstrip('.')
        return s


    major_locator = ticker.LogLocator(base=10., subs=subs)
    major_formatter = ticker.FuncFormatter(myfuncformatter)

    minor_locator = ticker.LogLocator(base=10., subs=[3., 4., 6., 7., 8., 9.])
    minor_formatter = ticker.NullFormatter()

    if 'x' in axis.lower():
        ax.xaxis.set_major_locator(major_locator)
        ax.xaxis.set_major_formatter(major_formatter)
        ax.xaxis.set_minor_locator(minor_locator)
        ax.xaxis.set_minor_formatter(minor_formatter)
    if 'y' in axis.lower():
        ax.yaxis.set_major_locator(major_locator)
        ax.yaxis.set_major_formatter(major_formatter)
        ax.yaxis.set_minor_locator(minor_locator)
        ax.yaxis.set_minor_formatter(minor_formatter)

    if grid:
        ax.grid(True, which="major", color=color, linestyle=":")
        ax.grid(True, which="minor", color=color, linestyle=":")


# ______________________________________________
def Ntick(ax, N, axis="xy"):
    if "x" in axis:
        ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=N))
    if "y" in axis:
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=N))


# ################################################ colors
def value2color(value, vmin=0., vmax=1.0, cmap = plt.cm.jet):
    index  = int(np.floor(cmap.N * (value - vmin) / (vmax - vmin)))
    return np.array(cmap(index)[:3])


# ______________________________________________
def values2colors(values, vmin=0., vmax=1.0, cmap = plt.cm.jet):
    indexs = np.floor(cmap.N * (values - vmin) / (vmax - vmin))
    indexs = np.array(indexs, int)
    return np.array([cmap(index)[:3] for index in indexs])


# ______________________________________________
def makecolorbar(vmin, vmax, cmap=plt.cm.get_cmap('jet')):
    cb = plt.cm.ScalarMappable(norm=None, cmap=cmap)
    cb.set_array([vmin, vmax])
    cb.set_clim((vmin, vmax))
    return cb

# ______________________________________________
def legendtext(ax, txt, fontsize = 14, loc = 2, multialignment = "left", frameon=True, framealpha=0.5, **kwargs):
    "create a legend with text only"
    at = AnchoredText(txt,
                      prop=dict(size=fontsize, multialignment = multialignment, **kwargs),
                      frameon=frameon,
                      loc=loc)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    at.patch.set_alpha(framealpha)
    ax.add_artist(at)
    return at


# _____________________________________
def textonly(ax, txt, fontsize = 14, loc = 2, multialignment = "left", frameon=True, framealpha=0.5, **kwargs):
    "create a legend with text only"
    at = AnchoredText(txt,
                      prop=dict(size=fontsize, multialignment = multialignment, **kwargs),
                      frameon=frameon,
                      loc=loc)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    at.patch.set_alpha(framealpha)
    ax.add_artist(at)
    return at

