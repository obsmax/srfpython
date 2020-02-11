from srfpython.standalone.multipro8 import Job, MapSync
from srfpython.standalone.stdout import waitbar
from srfpython.Herrmann.Herrmann import HerrmannCallerFromLists
from parameterizers import Parameterizer
from datacoders import Datacoder
import numpy as np

"""
Parameterizer : object to convert an array of parameters into smth that can be understood by Herrmann.dispersion
Theory        : object to convert a model array into a data array
Datacoder     : object to convert output from Herrmann.dispersion to an array of data

                          a model array (m)                        
  |                          ^    |
  |                          |    |
  |      mod96file ----->  parameterizer   --------> m_apr, CM
  |      (apriori)           |    |
  |                          |    v
  |                        depth model 
  |                     (ztop, vp, vs, rh)
  |                            |
 theory                 Herrmann.dispersion
 (forward problem)             |
  |                            v
  |                       dispersion data 
  |                      (waves, types, modes, freqs, values, (dvalues))
  |                          ^    |
  |                          |    |
  |     surf96file ----->   datacoder      --------> d_obs, CD
  |      (target)            |    |
  v                          |    v
                          a data array (d)
"""


class Theory(object):

    def __init__(self, parameterizer, datacoder, h=0.005, ddc=0.005):
        assert isinstance(parameterizer, Parameterizer)
        assert isinstance(datacoder, Datacoder)
        self.parameterizer, self.datacoder = parameterizer, datacoder

        self.herrmanncaller = HerrmannCallerFromLists(
            waves=datacoder.waves, types=datacoder.types,
            modes=datacoder.modes, freqs=datacoder.freqs,
            h=h, ddc=ddc)

    def __call__(self, m):
        ztop, vp, vs, rh = self.parameterizer.inv(m)   # recover model from parameterized array (m)

        values = self.herrmanncaller.disperse(ztop, vp, vs, rh)

        return self.datacoder(values)  # convert dispersion data to coded array  (d)


def overdisp(ms, overwaves, overtypes, overmodes, overfreqs, verbose=True, **mapkwargs):
    """extrapolate dispersion curves"""

    herrmanncaller = HerrmannCallerFromLists(
        waves=overwaves, types=overtypes,
        modes=overmodes, freqs=overfreqs,
        h=0.005, ddc=0.005)

    def fun(mms):
        ztop, vp, vs, rh = mms
        try:
            overvalues = herrmanncaller.dispersion(ztop, vp, vs, rh)

        except KeyboardInterrupt:
            raise

        except Exception:
            h = ztop[1:] - ztop[:-1]
            # assume failuer was caused by rounding issues
            h[h <= 0.001] = 0.001001
            ztop = np.concatenate(([0.], h.cumsum()))
            try:  # again
                overvalues = herrmanncaller.disperse(ztop, vp, vs, rh)

            except KeyboardInterrupt:
                raise

            except Exception:
                raise
                overvalues = np.nan * np.ones(len(overwaves))

        return mms, overvalues

    with MapSync(fun, (Job(mms) for mms in ms), **mapkwargs) as ma:
        if verbose: wb = waitbar('overdisp')
        Njobs = len(ms) - 1.
        for jobid, (mms, overvalues), _, _ in ma:
            if verbose: wb.refresh(jobid / Njobs)
            dds = (overwaves, overtypes, overmodes, overfreqs, overvalues)
            yield mms, dds
        if verbose:
            wb.close()
            print
    print
