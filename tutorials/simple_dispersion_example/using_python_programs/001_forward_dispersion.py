# -----------------------
# import all components of srfpython
# -----------------------
from srfpython import *

# -----------------------
# load a 1-D depth model created by 000_create_dephmodel.py
# -----------------------
dm = depthmodel_from_mod96('./model000.mod96')


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

results = dispersion_2(
    ztop=dm.vp.z, 
    vp=dm.vp.values, 
    vs=dm.vs.values, 
    rh=dm.rh.values, 
    Curves=Curves)  # just a generator no results generated until looping over it

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
