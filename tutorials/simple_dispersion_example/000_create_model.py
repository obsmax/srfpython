from srfpython import *
from srfpython.depthdisp.depthmodels import gardner74

# -----------------------
# create 1-D depth model using 4 arrays with same length
# -----------------------
if False:
    ztop = [0.00, 0.08, 0.45, 0.6, 1.1, 1.55, 2.33, 2.80]    # km, top layer depth, positive, increasing downward, 0 = surface
    vs   = [0.86, 1.30, 1.34, 1.69, 1.60, 2.53, 3.13, 3.31]  # km/s, S wave velocity in each layer
    pr   = [2.0,  1.9,  1.8,  1.85, 1.8,  1.79, 1.75, 1.73]  # VP/VS in each layer
    rh   = [2.27, 2.37, 2.47, 2.42, 2.52, 2.58, 2.58, 2.63]  # g/cm3, Density in each layer
else:
    # ztop = depthspace(2.8, 50)
    ztop = np.linspace(0., 2.8, 50)  # np.sort(np.random.rand(30) * 2.8)

    pr = 1.73 + .4 * np.exp(-ztop / 1.5) + \
         .01 * np.sin(2. * np.pi * ztop / 1.) + \
         .02 * np.sin(2. * np.pi * ztop / 2.) + \
         .005 * np.cos(2. * np.pi * ztop / 3.)

    vs = (3.5 - .86) / (ztop[-1] - ztop[0]) * (ztop - 0.) + .86 + \
         -.7 * np.exp(-ztop / .1) + \
         .08 * np.cos(2. * np.pi * ztop / .5) + \
         .2 * np.sin(2. * np.pi * ztop / 1.) + \
         .1 * np.cos(2. * np.pi * ztop / 2.) + \
         .15 * np.cos(2. * np.pi * ztop / 3.)

    # gardner74(pr * vs) #
    rh = (2.65 - 2.39) / (ztop[-1] - ztop[0]) * (ztop - 0.) + 2.39 + \
         .03 * np.sin(2. * np.pi * ztop / 1.) + \
         .01 * np.sin(2. * np.pi * ztop / 2.) + \
         .02 * np.cos(2. * np.pi * ztop / 3.)

# create the depthmodel object, use a subclass that is to be intitiated with arrays
# see also depthmodel, depthmodel1D, depthmodel_from_mod96, ...
dm = depthmodel_from_arrays(ztop, np.array(pr) * np.array(vs), vs, rh)

# __str__ returns the file content at mod96 format, (see Herrmann CPS documentation)
print dm 


dm.write96('model000.mod96')

