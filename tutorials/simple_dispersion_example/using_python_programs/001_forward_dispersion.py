# -----------------------
# import all components of srfpython
# -----------------------
from srfpython import *

# -----------------------
# load a 1-D depth model created by 000_create_dephmodel.py
# -----------------------
dm = depthmodel_from_mod96('./model000.mod96')


# __str__ returns the file content at mod96 format, (see Herrmann CPS documentation)
print dm 

# -----------------------
# compute dispersion curves from the depthmodel above
# -----------------------

# define the dipsersion curves to compute

curves = [Curve(wave='R', type='U', mode=0, freqs=freqspace(0.2, 3.5, 35, "log")),
          Curve(wave='R', type='U', mode=1, freqs=freqspace(0.2, 3.5, 35, "log")),
          Curve(wave='R', type='C', mode=0, freqs=freqspace(0.2, 3.5, 35, "log")),
          Curve(wave='R', type='C', mode=1, freqs=freqspace(0.2, 3.5, 35, "log")),
          Curve(wave='L', type='U', mode=0, freqs=freqspace(0.2, 3.5, 35, "log")),
          Curve(wave='L', type='U', mode=1, freqs=freqspace(0.2, 3.5, 35, "log")),
          Curve(wave='L', type='C', mode=0, freqs=freqspace(0.2, 3.5, 35, "log")),
          Curve(wave='L', type='C', mode=1, freqs=freqspace(0.2, 3.5, 35, "log"))]

# compute dispersion curves and display
hc = HerrmannCaller(curves=curves)
curves_out = hc(
    ztop=dm.vp.z,
    vp=dm.vp.values, 
    vs=dm.vs.values, 
    rh=dm.rh.values, 
    keepnans=False)

ax = plt.gca()
for curve in curves_out:
    curve.plot(ax, "+-")


ax.set_xscale('log')
ax.set_yscale('log')
logtick(ax, "xy")
ax.set_title('figure 2 : Herrmann.py demo')

plt.legend()
plt.show()
