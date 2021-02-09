from srfpython import *
from srfpython.sensitivitykernels.sker17 import sker17_1

"""
python script to compute depth sensitivity kernels
"""


# ================== build a 1d velocity model
ztop = np.arange(0., 30., 1.)  # top layer array, sorted, in km, 1st at 0.0
vs = np.linspace(1.8, 4.0, len(ztop))  # Vs in each layer from top to half space, km/s
vp = 1.73 * vs   # Vp in each layer, km/s
rh = 2.67 * np.ones_like(vs)  # Density in each layer g/cm3

# ================== compute sensitivity kernels for one wave type (RC0)
norm = True      # !!! recommended for non uniform layer model (try irregular ztop to see the difference) !!!
generator = sker17_1(ztop, vp, vs, rh,
    Waves=["R"],         # R=rayleigh
    Types=["C"],         # C=phase velocity
    Modes=[0],           # 0=fundamental mode
    Freqs=[[0.1, 0.2]],  # frequencies at which to compute 1d kernels, in Hz
    norm=norm)

# only one item here
wave, type, mode, freqs, DLOGVADZ, DLOGVADLOGVS, DLOGVADLOGPR, DLOGVADLOGRH = \
    generator.next()

# ================== display
plt.figure()

# depth model
ax1 = plt.subplot(121)
dm = depthmodel_from_arrays(ztop, vp, vs, rh)
dm.show(ax1)
ax1.grid(True)
plt.legend()

# kernels
if norm:
    title = r'$ \frac{H}{H_i} \, \frac{d lnV_{{%s%s%d}}}{d lnVs} $' % (wave, type, mode)
else:
    title = r'$ \frac{d lnV_{{%s%s%d}}}{d lnVs} $' % (wave, type, mode)

ax2 = plt.subplot(122, sharey=ax1, title=title, xlabel="sensitivity")
ax2.plot(DLOGVADLOGVS[:, 0], ztop, label="%s%s%d@%.2fHz" % (wave, type, mode, freqs[0]))
ax2.plot(DLOGVADLOGVS[:, 1], ztop, label="%s%s%d@%.2fHz" % (wave, type, mode, freqs[1]))
ax2.grid(True)
plt.legend()
plt.ion()
plt.show()
raw_input('pause')





