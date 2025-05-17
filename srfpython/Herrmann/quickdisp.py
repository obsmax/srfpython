#!/usr/bin/env python

import sys, glob, os
import numpy as np
import matplotlib.pyplot as plt
from srfpython.depthdisp.depthmodels import depthmodel_from_arrays
from srfpython.Herrmann.Herrmann import HerrmannCaller, Curve


ztopvpvsrh = np.asarray(sys.argv[1:], float)
nlayer = len(ztopvpvsrh) // 4
ztop = ztopvpvsrh[:nlayer]  #km
vp = ztopvpvsrh[nlayer:2*nlayer]  #km/s
vs = ztopvpvsrh[2*nlayer:3*nlayer]  #km/s
rh = ztopvpvsrh[3*nlayer:]  #g/cm3

dm = depthmodel_from_arrays(
	z=ztop, vp=vp, vs=vs, rh=rh)

cmin_approx = 0.9 * min(vs)
cmax_approx = 0.9 * max(vs)

lambda_min_approx = ztop[1] / 3.
lambda_max_approx = ztop[-1] * 3.

freq_min = cmin_approx / lambda_max_approx
freq_max = cmax_approx / lambda_min_approx

assert freq_min < freq_max 
freqs = np.logspace(np.log10(freq_min), np.log10(freq_max), 30)

print(f"{ztop=}")
print(f"{vp=}")
print(f"{dm.vs.stairs()=}")
print(f"{dm.rh.stairs()=}")
print(f"{freqs=}")


curves = [
    Curve(wave="R", type="C", mode=0, freqs=freqs),
    Curve(wave="R", type="C", mode=1, freqs=freqs),
    Curve(wave="R", type="C", mode=2, freqs=freqs),
    Curve(wave="R", type="C", mode=3, freqs=freqs),
    ]

hc = HerrmannCaller(curves=curves, h=0.0005, ddc=0.0005)
curves = hc(ztop=dm.vp.ztop(), vp=dm.vp.values, vs=dm.vs.values, rh=dm.rh.values)


fig = plt.figure()

ax2 = fig.add_subplot(222)
for curve in curves:
    ax2.plot(
        curve.freqs, curve.values * 1000., 
        label="%s%s%d" % (curve.wave, curve.type, curve.mode))
ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Phase Velocity (m/s)')

ax2.legend()

ax1 = fig.add_subplot(121)
dm.show(ax=ax1, units="m/s/kg.m-3")
ax1.set_ylabel('Depth (m)')
ax1.set_xlabel('Value (m/s, kg/m^3)')
ax1.legend()
plt.show()

    


