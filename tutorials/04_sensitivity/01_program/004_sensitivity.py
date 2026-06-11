import numpy as np
import matplotlib.pyplot as plt
from srfpython import depthmodel_from_arrays, Curve
from srfpython.Herrmann.Herrmann import HerrmannCaller
from srfpython.sensitivitykernels.sker17 import sker17
from srfpython.standalone.cmaps import tomocmap1
from srfpython.standalone.display import logtick, textonly

"""
experimental
"""

# ============
# build depth model : !!UNITS!! (km, km/s, km/s, g/cm3)
dm = depthmodel_from_arrays(
    z   =np.array([0.,      0.5,    1.5,    5.0,     10.,        ]) * 0.001,
    vp  =np.array([    3800.,  1800.,   1500.,  2100.,   3200.,  ]) * 0.001,
    vs  =np.array([    2030.,   600.,    500.,   800.,   1500.,  ]) * 0.001,
    rh  =np.array([    2400.,  2000.,   2000.,  2000.,   2300.,  ]) * 0.001,
    )

# add a layer in the top half space before depth resampling
dm.add_layer(thickness=10. * 0.001)  # in km

# split the continuous layers into smaller layers
dm.split(thickness=0.25 * 0.001)  # in km

if False:
    # display (qc)
    dm.show(plt.gca(), ".-", units='m/s/kg.m-3',)
    plt.show()

# define the expected dispersion curves
freqs = np.linspace(5., 150., 70) # Hz
curves = [
    Curve(wave="R", type="C", mode=0, freqs=freqs),
    Curve(wave="R", type="C", mode=1, freqs=freqs),
    Curve(wave="R", type="C", mode=2, freqs=freqs),
    Curve(wave="R", type="C", mode=3, freqs=freqs),
    ]

if False:
    # forward problem
    hc = HerrmannCaller(curves, ddc=0.005, h=0.005)
    curves = hc.call_dm(dm=dm)  # overwrite with output disp curves

    # display (qc)
    for curve_out in curves:
        plt.plot(
            curve_out.freqs,
            curve_out.values * 1000.,
            label="%s%s%d" % (curve_out.wave, curve_out.type, curve_out.mode))

    plt.legend()
    plt.show()

# sensitivity kernels
sker_gen = sker17(
    dm=dm,
    curves=curves,
    dlogz=0.01, dlogvp=0.01, dlogvs=0.01, dlogrh=0.01,
    norm=True, relative=False,
    h=0.005, ddc=0.005)

# ===================
for wave, typ, mode, skernels in sker_gen:

    for skernel in skernels:
        fig = skernel.show(vmin=None, vmax=None, cmap=tomocmap1(w=0.01, W=0.2), units="m/s/kg.m-3")
        fig.suptitle(f"{skernel.parameter_name}, {skernel.curve.wave}{skernel.curve.type}{skernel.curve.mode}")

plt.show()

# ===================
