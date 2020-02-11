from srfpython.Herrmann.Herrmann import HerrmannCaller, Curve
import numpy as np

nlayer = 10
nfreq = 30
ztop = np.linspace(0., 3., nlayer)
vs = np.linspace(1., 3.5, nlayer)
vp = 1.73 * vs
rh = np.linspace(2.5, 3.1, nlayer)
f = np.logspace(-1, 1, nfreq)

curves = [Curve(wave="R", type="C", mode=0, freqs=f),
          Curve(wave="R", type="U", mode=0, freqs=f),
          Curve(wave="L", type="C", mode=0, freqs=f),
          Curve(wave="L", type="U", mode=0, freqs=f)]

hc = HerrmannCaller(curves=curves, h=0.005, ddc=0.005)

for _ in range(100):
    hc.disperse(ztop, vp, vs, rh)
