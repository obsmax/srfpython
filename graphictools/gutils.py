import matplotlib.pyplot as plt
from matplotlib import _pylab_helpers

gcf, gca = plt.gcf, plt.gca
# -------------------------
def pause(message = "pause"):
    raw_input(message)
# -------------------------
def gcfs():
    "get current figures"
    return [w.canvas.figure for w in _pylab_helpers.Gcf.get_all_fig_managers()]
# -------------------------
def showme(block = True):
    figures = gcfs()
    for figure in figures:
        figure.canvas.draw()
        figure.show()
    if block:
        pause()