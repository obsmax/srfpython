import numpy as np
from srfpython.Herrmann.Herrmann import HerrmannCallerFromLists
from srfpython.standalone.multipro8 import Job, MapSync
from srfpython.standalone.stdout import waitbar


class _OverdispCore(object):
    def __init__(self, herrmanncaller):
        self.herrmanncaller = herrmanncaller

    def __call__(self, mms):
        ztop, vp, vs, rh = mms
        try:
            overvalues = self.herrmanncaller.dispersion(ztop, vp, vs, rh)

        except KeyboardInterrupt:
            raise

        except Exception as e:
            # assume failure was caused by rounding issues

            h = ztop[1:] - ztop[:-1]
            h[h <= 0.001] = 0.001001
            ztop = np.concatenate(([0.], h.cumsum()))
            try:  # again
                overvalues = self.herrmanncaller.disperse(ztop, vp, vs, rh)

            except:
                raise e  # raise the initial error

        return mms, overvalues


def overdisp(ms, overwaves, overtypes, overmodes, overfreqs, verbose=True, **mapkwargs):
    """extrapolate dispersion curves"""

    herrmanncaller = HerrmannCallerFromLists(
        waves=overwaves, types=overtypes,
        modes=overmodes, freqs=overfreqs,
        h=0.005, ddc=0.005)

    fun = _OverdispCore(herrmanncaller)
    gen = (Job(mms) for mms in ms)

    with MapSync(fun, gen, **mapkwargs) as ma:
        if verbose: wb = waitbar('overdisp')
        Njobs = len(ms) - 1.
        for jobid, (mms, overvalues), _, _ in ma:
            if verbose: wb.refresh(jobid / Njobs)
            dds = (overwaves, overtypes, overmodes, overfreqs, overvalues)
            yield mms, dds
        if verbose:
            wb.close()
            print()
    print()