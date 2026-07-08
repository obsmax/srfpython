#!/usr/bin/env python

"""
Display Brocher 2005 empirical relations for vp(vs) and rho(vp(vs))
"""
import numpy as np
import matplotlib.pyplot as plt
from srfpython.depthdisp.depthmodels import brocher2005


vs = np.linspace(50., 1000., 1000)

vp, rh = brocher2005(vs/1000.)
vp *= 1000.
rh *= 1000.
vpavs = vp / vs
nu = 0.5 * (1 - 1. / ((vpavs)**2. - 1.))

fig = plt.figure()
fig.suptitle('Brocher 2005 empirical relations')
plt.subplot(311, xlabel="vs")
plt.plot(vs, vp, label="vp(vs)")
plt.plot(vs, rh, label="rho(vp(vs))")
plt.gca().legend()

plt.subplot(312, sharex=plt.gca(), xlabel="vs")
plt.plot(vs, vpavs, label="vp(vs)/vs")
plt.gca().legend()

plt.subplot(313, sharex=plt.gca(), xlabel="vs")
plt.plot(vs, nu, label=r"$\nu$")
plt.gca().legend()
plt.show()
#plt.plot(vs, vp/vs, label="vp/vs")