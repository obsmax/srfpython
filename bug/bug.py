from srfpython import *

for i in range(1):
    if i == 0:
        f = freqspace(0.2, 2.0, 20, "plog")
    elif i == 1:
        f = freqspace(0.2, 2.0, 200, "plog")
    curves = [
        Curve(wave="R", type="U", mode=0, freqs=f),
        Curve(wave="R", type="U", mode=1, freqs=f)]

    hc = HerrmannCaller(curves=curves, h=0.005, ddc=0.005)
    dm = depthmodel_from_mod96('node014.mod96')
    curves_out = hc(dm.vs.z, dm.vp.values, dm.vs.values, dm.rh.values, keepnans=False)

    for curve in curves_out:
        curve.plot(plt.gca(), 'o-')

plt.ion()
plt.show()
input('')
