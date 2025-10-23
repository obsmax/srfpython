#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt


vpavs = np.linspace(1.41424, 20., 10000)
nu = 0.5 * (1 - 1. / ((vpavs)**2. - 1.))

plt.plot(nu, vpavs)
plt.xlabel(r"$ \nu $")
plt.ylabel(r"$ \frac{V_p}{V_s} $")
plt.grid(True)
plt.show()
