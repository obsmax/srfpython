from __future__ import print_function
# -----------------------
# import all components of srfpython
# -----------------------
from srfpython import *

# -----------------------
# load a 1-D depth model created by 000_create_dephmodel.py
# -----------------------
dm = depthmodel_from_mod96('./model000.mod96')


# __str__ returns the file content at mod96 format, (see Herrmann CPS documentation)
print(dm)

# -----------------------
# compute dispersion curves from the depthmodel above
# -----------------------

# define the dipsersion curves to compute
f = freqspace(0.2, 3.5, 35, "log")
curves = [Curve(wave='R', type='U', mode=0, freqs=f),
          Curve(wave='R', type='U', mode=1, freqs=f),
          Curve(wave='R', type='C', mode=0, freqs=f),
          Curve(wave='R', type='C', mode=1, freqs=f),
          Curve(wave='L', type='U', mode=0, freqs=f),
          Curve(wave='L', type='U', mode=1, freqs=f),
          Curve(wave='L', type='C', mode=0, freqs=f),
          Curve(wave='L', type='C', mode=1, freqs=f)]

# compute dispersion curves and display
hc = HerrmannCaller(curves=curves)
curves_out = hc(
    ztop=dm.vp.z,
    vp=dm.vp.values, 
    vs=dm.vs.values, 
    rh=dm.rh.values, 
    keepnans=False)

fig = plt.figure()
ax = fig.add_subplot(111, xscale="log", yscale="log")
for curve in curves_out:
    curve.plot(ax, "+-")

logtick(ax, "xy")
ax.set_title('figure 2 : Herrmann.py demo')

plt.legend()
plt.show()
