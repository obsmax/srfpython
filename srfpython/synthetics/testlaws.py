import numpy as np
from srfpython.depthdisp.dispcurves import mklaws, igroupbywtm
from srfpython.depthdisp.depthmodels import depthmodel_from_arrays
from srfpython.Herrmann.Herrmann import HerrmannCaller, Curve

# depth model
ztop = [0.00, 0.25, 0.45, 0.65, 0.85, 1.05, 1.53, 1.80]  # km, top layer depth
vp = [1.85, 2.36, 2.63, 3.15, 3.71, 4.54, 5.48, 5.80]    # km/s
vs = [0.86, 1.10, 1.24, 1.47, 1.73, 2.13, 3.13, 3.31]    # km/s
rh = [2.47, 2.47, 2.47, 2.47, 2.47, 2.58, 2.58, 2.63]    # g/cm3

dm = depthmodel_from_arrays(ztop, vp, vs, rh)

# dipsersion parameters
# f = np.logspace(np.log10(0.1), np.log10(500.), 100) # frequency array, km/s : ERROR : example of missing modes as reported by Herrmann !!!!!
f = np.logspace(np.log10(0.1), np.log10(10.), 100) # frequency array, km/s

# dispersion curves
curves = [Curve(wave='R', type='C', mode=0, freqs=f),
          Curve(wave='R', type='C', mode=1, freqs=f),
          Curve(wave='R', type='C', mode=2, freqs=f),
          Curve(wave='R', type='C', mode=3, freqs=f),
          Curve(wave='R', type='C', mode=4, freqs=f),
          Curve(wave='R', type='C', mode=5, freqs=f),
          Curve(wave='R', type='C', mode=6, freqs=f)]

# compute dispersion laws
#Waves, Types, Modes, Freqs = zip(*curves)
#waves, types, modes, freqs = igroupbywtm(Waves, Types, Modes, Freqs)
hc = HerrmannCaller(curves)
#values = dispersion(ztop, vp, vs, rh, waves, types, modes, freqs)
values = hc.disperse(ztop, vp, vs, rh)
laws = mklaws(hc.waves, hc.types, hc.modes, hc.freqs, values, dvalues=None)
c0, c1, c2, c3 = laws[:4]

if __name__ == "__main__":
    from srfpython import *

    plt.figure()
    dm.show(gca())
    plt.legend()

    plt.figure()
    for law in laws:
        law.show(gca(), period=True)
    plt.legend()
    showme()