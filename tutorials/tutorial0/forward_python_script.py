# -----------------------
# import all components of srfpython
# -----------------------
from srfpython import *

# -----------------------
# create 1-D depth model using 4 arrays with same length
# -----------------------
ztop = [0.00, 0.25, 0.45, 0.65, 0.85, 1.05, 1.53, 1.80]  # km, top layer depth, positive, increasing downward, 0 = surface
vp   = [1.85, 2.36, 2.63, 3.15, 3.71, 4.54, 5.48, 5.80]  # km/s, P wave velocity in each layer
vs   = [0.86, 1.10, 1.24, 1.47, 1.73, 2.13, 3.13, 3.31]  # km/s, S wave velocity in each layer
rh   = [2.47, 2.47, 2.47, 2.47, 2.47, 2.58, 2.58, 2.63]  # g/cm3, Density in each layer

# create the depthmodel object, use a subclass that is to be intitiated with arrays
# see also depthmodel, depthmodel1D, depthmodel_from_mod96, ...
dm = depthmodel_from_arrays(ztop, vp, vs, rh)

# __str__ returns the file content at mod96 format, (see Herrmann CPS documentation)
print dm 


# -----------------------
# create 1-D depth model using 4 arrays with same length
# -----------------------
ztop = [0.00, 0.25, 0.45, 0.65, 0.85, 1.05, 1.53, 1.80]  # km, top layer depth, positive, increasing downward, 0 = surface
vp   = [1.85, 2.36, 2.63, 3.15, 3.71, 4.54, 5.48, 5.80]  # km/s, P wave velocity in each layer
vs   = [0.86, 1.10, 1.24, 1.47, 1.73, 2.13, 3.13, 3.31]  # km/s, S wave velocity in each layer
rh   = [2.47, 2.47, 2.47, 2.47, 2.47, 2.58, 2.58, 2.63]  # g/cm3, Density in each layer

# create the depthmodel object, use a subclass that is to be intitiated with arrays
# see also depthmodel, depthmodel1D, depthmodel_from_mod96, ...
dm = depthmodel_from_arrays(ztop, vp, vs, rh)


# __str__ returns the file content at mod96 format, (see Herrmann CPS documentation)
print dm 

# -----------------------
# compute dispersion curves from the depthmodel above
# -----------------------

# define the dipsersion curves to compute
#          Wave(R/L) Type(C/U) Mode    Frequency array (Hz)             
Curves = [('R',      'U',      0,      freqspace(0.2, 3.5, 35, "log")), 
          ('R',      'U',      1,      freqspace(0.2, 3.5, 35, "log")), 
          ('R',      'C',      0,      freqspace(0.2, 3.5, 35, "log")), 
          ('R',      'C',      1,      freqspace(0.2, 3.5, 35, "log")), 
          ('L',      'U',      0,      freqspace(0.2, 3.5, 35, "log")), 
          ('L',      'U',      1,      freqspace(0.2, 3.5, 35, "log")), 
          ('L',      'C',      0,      freqspace(0.2, 3.5, 35, "log")), 
          ('L',      'C',      1,      freqspace(0.2, 3.5, 35, "log"))] 

# compute dispersion curves and display
results = dispersion_2(ztop, vp, vs, rh, Curves)  # just a generator no results generated until looping over it

ax = plt.gca()
for wave_letter, type_letter, mode_number, \
    frequency_array, velocity_array in results:

    ax.loglog(1. / frequency_array, velocity_array, '+-',
              label="%s%s%d" % (wave_letter, type_letter, mode_number))

ax.set_xlabel('period (s)')
ax.set_ylabel('velocity (km/s)')    
ax.grid(True, which="major")
ax.grid(True, which="minor")
logtick(ax, "xy")
ax.set_title('figure 2 : Herrmann.py demo')

plt.legend()
plt.show()
