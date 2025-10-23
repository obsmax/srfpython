import numpy as np
import matplotlib.pyplot as plt

from srfpython.HerrMet.plugins.optimize import HERRMETOPTIMIZEPRIORFILE, HERRMETOPTIMIZEDOBSFILE, HERRMETOPTIMIZEINVFILE
from srfpython.HerrMet.paramfile import load_paramfile
from srfpython.HerrMet.datacoders import makedatacoder, Datacoder_log
from srfpython.depthdisp.depthmodels import depthmodel_from_arrays
from srfpython.depthdisp.depthdispdisplay import DepthDispDisplay


with np.load(HERRMETOPTIMIZEPRIORFILE) as loader:
    mprior = loader['Mprior'][:, 0, 0]
    parameterizer_string = loader['parameterizer_strings'][0]
    parameterizer, logRHOM = load_paramfile(parameterizer_string, verbose=True)

with np.load(HERRMETOPTIMIZEDOBSFILE) as loader:
    datacoder_string = loader['datacoder_strings'][0]
    datacoder = makedatacoder(datacoder_string, which=Datacoder_log)
    dobs = loader['Dobs']

with np.load(HERRMETOPTIMIZEINVFILE) as loader:
    m0 = loader['M0'][:, 0, 0]
    msol = loader['Msol'][:, 0, 0]
    d0 = loader['D0']
    dsol = loader['Dsol']

dm_prior = depthmodel_from_arrays(*parameterizer.inv(mprior))
dm_0 = depthmodel_from_arrays(*parameterizer.inv(m0))
dm_sol = depthmodel_from_arrays(*parameterizer.inv(msol))

fig = plt.figure()
ddd = DepthDispDisplay(fig=fig, targetfile=datacoder_string)

ddd.plotdisp(waves=datacoder.waves, types=datacoder.types, modes=datacoder.modes, freqs=datacoder.freqs,
             values=datacoder.inv(dobs),
             color="m", linewidth=3, alpha=1, label="obs")

ddd.plotdisp(waves=datacoder.waves, types=datacoder.types, modes=datacoder.modes, freqs=datacoder.freqs,
             values=datacoder.inv(d0),
             color="k", linewidth=1, alpha=1, label="start")

ddd.plotdisp(waves=datacoder.waves, types=datacoder.types, modes=datacoder.modes, freqs=datacoder.freqs,
             values=datacoder.inv(dsol),
             color="r", linewidth=3, alpha=1, label="optimized")


ddd.plot_depthmodel(dm_prior, color="m", linewidth=3, alpha=1, label="prior")
ddd.plot_depthmodel(dm_0, color="k", linewidth=1, alpha=1, label="start")
ddd.plot_depthmodel(dm_sol, color="r", linewidth=3, alpha=1, label="optimized")
fig.axes[0].legend()
plt.show()
