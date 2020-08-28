import numpy as np
import matplotlib.pyplot as plt
from srfpython.depthdisp.dispcurves import surf96reader
from srfpython.HerrMet.datacoders import Datacoder_log, makedatacoder
from srfpython.standalone.stdout import waitbarpipe

import scipy.sparse as sp

# ==================== read the parameter file
with open('optimize.param', 'r') as fid:
    for l in fid:
        l = l.split('\n')[0].split('#')[0].strip()
        if not len(l):
            continue
        if l.startswith('nx'):           nx = int(l.split()[1])
        if l.startswith('ny'):           ny = int(l.split()[1])
        if l.startswith('newnz'):     newnz = int(l.split()[1])
        if l.startswith('ndecim'):   ndecim = int(l.split()[1])

        if l.startswith('dx'):           dx = float(l.split()[1])
        if l.startswith('dy'):           dy = float(l.split()[1])
        if l.startswith('newdz'):     newdz = float(l.split()[1])
        if l.startswith('Lh'):           Lh = float(l.split()[1])
        if l.startswith('Lv'):           Lv = float(l.split()[1])
        if l.startswith('vsunc'):     vsunc = float(l.split()[1])
        if l.startswith('damping'): damping = float(l.split()[1])

        if l.startswith('extractfiles'): extractfiles = l.split()[1]
        if l.startswith('datafiles'):       datafiles = l.split()[1]


x = (np.arange(nx) * dx)[::ndecim]
y = (np.arange(ny) * dy)[::ndecim]
ixs = np.arange(nx)[::ndecim]  # indexs used to name the files
iys = np.arange(ny)[::ndecim]


def get_datacoders(verbose=True):
    datacoder_strings = []

    for i, iy in enumerate(iys):
        for j, ix in enumerate(ixs):  # order matters!!!!

            try:
                datafile = datafiles.format(iy=iy, ix=ix)
                print('loading ', datafile)
                s96 = surf96reader(filename=datafile)


            except Exception as e:
                raise e

            datacoder_string = str(s96)
            datacoder_strings.append(datacoder_string)

    datacoder_strings = np.asarray(datacoder_strings, str)
    datacoders = [makedatacoder(datacoder_string, which=Datacoder_log) for datacoder_string in datacoder_strings]
    return datacoders, datacoder_strings


datacoders, datacoder_strings = get_datacoders()
Nobs = []
Waves = []
Types = []
Modes = []
Freqs = []
Dobs = []
Dunc = []
for datacoder in datacoders:
    _dobs, _CDinv = datacoder.target()
    Nobs.append(len(datacoder.waves))
    Waves.append(datacoder.waves)
    Types.append(datacoder.types)
    Modes.append(datacoder.modes)
    Freqs.append(datacoder.freqs)
    Dobs.append(_dobs)
    Dunc.append(_CDinv ** -0.5)

Nobs = np.array(Nobs)
Waves = np.hstack(Waves)
Types = np.hstack(Types)
Modes = np.hstack(Modes)
Freqs = np.hstack(Freqs)
Dobs = np.hstack(Dobs)
Dunc = np.hstack(Dunc)

np.savez('dobs.npz',
    datacoder_strings=datacoder_strings,
    Nobs=Nobs,
    Waves=Waves,
    Types=Types,
    Modes=Modes,
    Freqs=Freqs,
    Dobs=Dobs,
    Dunc=Dunc)

plt.figure()
plt.plot(Dobs)
plt.show()

# assert len(Dobs) == nx * ny * nper
# # CD = sp.diags(Dunc ** 2.0, shape=(len(Dobs), len(Dobs)), format="csc")
# # sp.save_npz('CD.npz', CD)
#
# np.save('datacoder_strings.npy', datacoder_strings)
# np.save('npernynx.npy', np.array([nper, ny, nx]))  # warning Dobs in a flat version of (ny, nx, nper)
# np.save('periods.npy', periods)
# np.save('pflat.npy', pflat)
# np.save('Dobs.npy', Dobs)
# np.save('Dunc.npy', Dunc)
