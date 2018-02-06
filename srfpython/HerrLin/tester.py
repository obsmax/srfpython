from srfpython import *

for i in xrange(10):
    ztop = np.sort(np.unique(np.random.rand(10) * 3.))
    ztop[0] = 0.

    vp = np.linspace(1., 5., len(ztop))
    pr = np.linspace(2.15, 1.8, len(ztop))
    rh = np.linspace(2., 3., len(ztop))
    vs = vp / pr

    dm = depthmodel_from_arrays(ztop, vp, vs, rh)
    dm.show(gca())
    dm.write96('%03d.mod96' % i)

    Curves = [('R', 'U', 0, freqspace(.2, 1., 15, "log")), \
              ('R', 'U', 1, freqspace(.2, 1., 15, "log")), \
              ('L', 'U', 0, freqspace(.2, 1., 15, "log")), \
              ('L', 'U', 1, freqspace(.2, 1., 15, "log"))]

    waves, types, modes, freqs, values = [[] for _ in xrange(5)]
    for w, t, m, F, V in dispersion_2(ztop, vp, vs, rh, Curves,
        h=0.005, dcl=0.005, dcr=0.005, keepnans=False):

        waves = np.concatenate((waves, np.array([w]).repeat(len(F))))
        types = np.concatenate((types, np.array([t]).repeat(len(F))))
        modes = np.concatenate((modes, np.array([m]).repeat(len(F))))
        freqs = np.concatenate((freqs, F))
        values = np.concatenate((values, V))

    s = surf96reader_from_arrays(waves, types, modes, freqs, values)
    s.write96('%03d.surf96' % i)

