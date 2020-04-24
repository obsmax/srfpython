from __future__ import print_function
from srfpython import *

# -----------------------
# create 1-D depth model using 4 arrays with same length
# -----------------------
# top layers array first at 0, positive, growing, km
ztop = np.linspace(0., 2.8, 50)

# vs in km/s
vs = (3.5 - .86) / (ztop[-1] - ztop[0]) * (ztop - 0.) + .86 + \
     -.7 * np.exp(-ztop / .1) + \
     .08 * np.cos(2. * np.pi * ztop / .5) + \
     .2 * np.sin(2. * np.pi * ztop / 1.) + \
     .1 * np.cos(2. * np.pi * ztop / 2.) + \
     .15 * np.cos(2. * np.pi * ztop / 3.)

vp, rh = brocher2005(vs)

# create the depthmodel object, use a subclass that is to be intitiated with arrays
# see also depthmodel, depthmodel1D, depthmodel_from_mod96, ...
dm = depthmodel_from_arrays(ztop, vp, vs, rh)

# __str__ returns the file content at mod96 format, (see Herrmann CPS documentation)
print(dm)

dm.write96('model000.mod96')
