import numpy as np
from srfpython.depthdisp.depthmodels import depthmodel_from_arrays

dm = depthmodel_from_arrays(
    z=np.array([0., 1., 30., 50.]),
    vp=np.array([1000., 1200., 1300., 1500.]),
    vs=np.array([ 800.,  900., 1000., 1100.]),
    rh=np.array([2200., 2300., 2500., 2600.]),
    qp=np.array([100., 200., 300., 400.]),
    qs=np.array([100., 200., 300., 400.]),    
    etap=np.array([0., 0., 0., 0.]),
    etas=np.array([0., 0., 0., 0.]),
    fp=np.array([1., 1., 1., 1.]),
    fs=np.array([1., 1., 1., 1.]),
    )
    
print(dm)
    
