#!/usr/bin/env bash

rm -f ./_HerrMet_data010/*rank*mod96
HerrMet --extract _HerrMet_data010 -top 3 0 1

#m96 --show model000.mod96 _HerrMet_data010/*rank*.mod96

/usr/bin/env python << END
import glob
from srfpython import *
dms = []
for f in glob.iglob('./_HerrMet_data010/*rank*mod96'):
    dm = depthmodel_from_mod96(f)
    dms.append(dm)
ztop = np.unique(np.concatenate([dm.ztop() for dm in dms]))
if np.any(ztop[1:] - ztop[:-1] < 0.001):
    ztop = np.arange(ztop[0], ztop[-1], 0.05)
dms = [dm.interp(ztop) for dm in dms]

vsmean = np.mean([dm.vs.values for dm in dms], axis=0)
prmean = np.mean([dm.pr().values for dm in dms], axis=0)
rhmean = np.mean([dm.rh.values for dm in dms], axis=0)
vpmean = prmean * vsmean
dmmean = depthmodel_from_arrays(ztop, vpmean, vsmean, rhmean)
dmmean.write96('_HerrMet_data010/mean.mod96')

END

m96 --disp  _HerrMet_data010/*{rank,mean}*.mod96 -RC0 .2 2. 100 plog -RU0 .2 2. 100  plog
