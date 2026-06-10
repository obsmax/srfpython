import numpy as np
import matplotlib.pyplot as plt

from srfpython.depthdisp.depthmodels import depthmodel_from_arrays
from srfpython.sensitivitykernels.sker17 import sker17, Curve, tomocmap1

"""
python script to compute depth sensitivity kernels
"""


# ================== build a 1d velocity model
dm = depthmodel_from_arrays(
    z = np.arange(0., 30., 1.),  # top layer array, sorted, in km, 1st at 0.0
    vp = 1.73 * np.linspace(1.8, 4.0, 30),   # Vp in each layer, km/s
    vs = np.linspace(1.8, 4.0, 30),  # Vs in each layer from top to half space, km/s
    rh = 2.67 * np.ones(30),  # Density in each layer g/cm3
    )

# ================== compute sensitivity kernels for one wave type (RC0)
curves = [Curve(wave='R', type='C', mode=0, freqs=np.linspace(0.005, 0.2, 50)),
          ]

sker_gen = sker17(
    dm=dm,
    curves=curves,
    norm=True, relative=False,
    ddc=0.005,
    )

# only one item here
# ===================
for wave, typ, mode, skernels in sker_gen:

    for skernel in skernels:
        fig = skernel.show(vmin=None, vmax=None, cmap=tomocmap1(w=0.01, W=0.2), units="km/s/g.cm-3")
        fig.suptitle(f"{skernel.parameter_name}, {skernel.curve.wave}{skernel.curve.type}{skernel.curve.mode}")
    plt.show()
