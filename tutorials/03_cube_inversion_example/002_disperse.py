import sys, glob, os
from srfpython.depthdisp.depthmodels import depthmodel_from_mod96
from srfpython.depthdisp.dispcurves import surf96reader_from_curves
from srfpython.Herrmann.Herrmann import Curve, HerrmannCaller
import numpy as np

os.system('mkdir --parents ./data')


freqs = np.logspace(np.log10(0.1), np.log10(1.), 10)


hc = HerrmannCaller(curves=[
        Curve(wave="R", type="U", mode=0, freqs=freqs),
        Curve(wave="R", type="C", mode=0, freqs=freqs),
        Curve(wave="R", type="U", mode=1, freqs=freqs),
        Curve(wave="R", type="C", mode=1, freqs=freqs)])

for fin in glob.glob('./model/nodes/*mod96'):
    dm = depthmodel_from_mod96(fin)
    curves = hc(ztop=dm.vs.z, vp=dm.vp.values, vs=dm.vs.values, rh=dm.rh.values, keepnans=False)
    s96r = surf96reader_from_curves(curves)

    fout = fin.replace('./model/nodes/', './data/').replace('.mod96', '.surf96')
    print(fin, '->', fout)
    s96r.write96(filename=fout)
