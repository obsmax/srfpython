from srfpython.Herrmann.Herrmann import dispersion_1
import numpy as np

nlayer = 10
nfreq = 30
ztop = np.linspace(0., 3., nlayer)
vs = np.linspace(1., 3.5, nlayer)
vp = 1.73 * vs
rh = np.linspace(2.5, 3.1, nlayer)
f = np.logspace(-1, 1, nfreq)
Waves = ['R', 'R', 'L', 'L']
Types = ['C', 'U', 'C', 'U']
Modes = [ 0,   0,   0,   0]
Freqs = [ f,   f,   f,   f]

for _ in range(100):

    for w, t, m, F, V in dispersion_1(ztop, vp, vs, rh, 
        Waves, Types, Modes, Freqs,
        h = 0.005, dcl = 0.005, 
        dcr = 0.005, keepnans = False):
        
        pass
        # print(w, t, m)
        
        
