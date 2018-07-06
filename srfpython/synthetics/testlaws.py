import numpy as np
from srfpython.depthdisp.dispcurves import mklaws, igroupbywtm
from srfpython.depthdisp.depthmodels import depthmodel_from_arrays
from srfpython.Herrmann.Herrmann import dispersion

# depth model
ztop = [0.00, 0.25, 0.45, 0.65, 0.85, 1.05, 1.53, 1.80]  # km, top layer depth
vp = [1.85, 2.36, 2.63, 3.15, 3.71, 4.54, 5.48, 5.80]  # km/s
vs = [0.86, 1.10, 1.24, 1.47, 1.73, 2.13, 3.13, 3.31]  # km/s
rh = [2.47, 2.47, 2.47, 2.47, 2.47, 2.58, 2.58, 2.63]  # g/cm3
dm = depthmodel_from_arrays(ztop, vp, vs, rh)

# dipsersion parameters
f = np.logspace(np.log10(0.2), np.log10(3.5), 35) # frequency array, km/s

# dispersion curves
curves = [('R', 'C', 0, f),
          ('R', 'C', 1, f),
          ('R', 'C', 2, f),
          ('R', 'C', 3, f),
          ('R', 'C', 4, f),
          ('R', 'C', 5, f),
          ('R', 'C', 6, f)]

# compute dispersion laws
Waves, Types, Modes, Freqs = zip(*curves)
waves, types, modes, freqs = igroupbywtm(Waves, Types, Modes, Freqs)
values = dispersion(ztop, vp, vs, rh, waves, types, modes, freqs)
laws = mklaws(waves, types, modes, freqs, values, dvalues=None)
c0, c1, c2, c3 = laws[:4]
