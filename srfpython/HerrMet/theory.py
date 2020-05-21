import warnings
import numpy as np
from srfpython.Herrmann.Herrmann import HerrmannCallerBasis
from srfpython.HerrMet.parameterizers import Parameterizer
from srfpython.HerrMet.datacoders import Datacoder

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
 theory                  Herrmann.HerrmannCaller.disperse
 (forward problem)             |
 (Frechet derivatives)         v
  |                       dispersion data 
  |                      (waves, types, modes, freqs, values, (dvalues))
  |                          ^    |
  |                          |    |
  |     surf96file ----->   datacoder      --------> d_obs, CD => logRHOD
  |      (target)            |    |
  v                          |    v
                          a data array (d)
"""


class Theory(object):

    def __init__(self, parameterizer, datacoder, h=0.005, ddc=0.005):
        if not isinstance(parameterizer, Parameterizer):
            raise TypeError(type(parameterizer))
        if not isinstance(datacoder, Datacoder):
            raise TypeError(type(datacoder))
        self.parameterizer, self.datacoder = parameterizer, datacoder

        self.herrmanncaller = HerrmannCallerBasis(
            waves=datacoder.waves, types=datacoder.types,
            modes=datacoder.modes, freqs=datacoder.freqs,
            h=h, ddc=ddc)

        # the output of self.herrmanncaller.disperse must match excatcly the datacoder order
        assert self.herrmanncaller.waves is datacoder.waves
        assert self.herrmanncaller.types is datacoder.types
        assert self.herrmanncaller.modes is datacoder.modes
        assert self.herrmanncaller.freqs is datacoder.freqs

    def __call__(self, m):
        """solves the forward problem"""

        # recover model from parameterized array (m)
        ztop, vp, vs, rh = self.parameterizer.inv(m)

        # call Herrmann's codes
        values = self.herrmanncaller.disperse(ztop, vp, vs, rh)

        # convert and return dispersion data to coded array (d)
        return self.datacoder(values)

    def frechet_derivatives(self, m, gm=None):
        """
        :param m: the model near which to compute frechet derivatives
        :param gm: if the dispersion data is known already (otherwise, I compute it using self.__call__(m)
        :return fd: the frechet derivatives
        """

        # compute the data for the current model if not provided
        if gm is None:
            gm = self.__call__(m)

        if np.isnan(gm).any():
            # probably because some freq points are above the cut off period
            warnings.warn('g(m) contains nans')
            # raise ValueError('g(m) contains nans')

        if np.isinf(gm).any():
            raise ValueError('g(m) contains infs')

        # shift the current model parameters by a delta
        deltam = self.parameterizer.frechet_deltas()
        if not np.all(deltam > 0.):
            raise ValueError('deltam must be positive')

        fd = np.zeros((len(gm), len(m)), float)
        for j in range(len(m)):
            mj = m.copy()
            mj[j] += deltam[j]
            gmj = self.__call__(mj)
            if np.isnan(gmj).any() or np.isinf(gmj).any():
                # probably because some freq points are above the cut off period
                warnings.warn('g(m+dm) contains nans of infs')
                # raise ValueError('g(m+dm) contains nans of infs')

            fd[:, j] = (gmj - gm) / deltam[j]

        fd[np.isnan(fd)] = 0.
        return fd


if __name__ == '__main__':
    from srfpython.standalone.asciifile import AsciiFile_fromstring
    from srfpython.HerrMet.parameterizers import Parameterizer_mZVSVPvsRHvp
    from srfpython.HerrMet.datacoders import Datacoder_log, makedatacoder

    if False:
        parameter_file_string = """
#met  NLAYER = 9
#met  TYPE = 'mZVSPRRH'
#met  PRIORTYPE = "DVPDVSDRHDPR"
#met DVSMIN = -0.5
#met DVSMAX = +1.5
#met DVPMIN = -0.5
#met DVPMAX = +1.5
#met DPRMIN = -1.0
#met DPRMAX = +0.0
#met DRHMIN = +0.0
#met DRHMAX = +1.0

#fld KEY     VINF          VSUP
#unt []      []            []
#fmt %s      %f            %f
     -Z1     -0.200000     -0.05
     -Z2     -0.400000     -0.2
     -Z3     -0.800000     -0.400000
     -Z4     -1.000000     -0.800000
     -Z5     -1.500000     -1.000000
     -Z6     -2.000000     -1.500000
     -Z7     -2.500000     -2.000000
     -Z8     -3.000000     -2.500000

     VS0     0.2           1.500000
     VS1     0.2           1.700000
     VS2     0.5           2.000000
     VS3     0.5           2.500000
     VS4     1.0           3.000000
     VS5     1.0           3.250000
     VS6     1.3           3.500000
     VS7     1.3           3.550000
     VS8     3.2           3.600000

     PR0     1.800000      3.6
     PR1     1.800000      2.3
     PR2     1.800000      2.3
     PR3     1.800000      2.3
     PR4     1.700000      2.3
     PR5     1.700000      2.3
     PR6     1.600000      2.0
     PR7     1.600000      2.0
     PR8     1.600000      1.9

     RH0     1.80          2.3
     RH1     1.80          2.3
     RH2     1.90          2.4
     RH3     2.00          2.4
     RH4     2.00          2.5
     RH5     2.10          2.5
     RH6     2.30          2.7
     RH7     2.40          2.7
     RH8     2.50          2.7
    """
    else:
        parameter_file_string = """
# BROCHER2005
#met RHvp = 'def RHvp(VP):return  1.6612 * VP - 0.4721 * VP ** 2.0 + 0.0671 * VP ** 3.0 - 0.0043 * VP ** 4.0 + 0.000106 * VP ** 5.0'
#met NLAYER = 10
#met TYPE = 'mZVSVPvsRHvp'
#met VPvs = 'def VPvs(VS):return  0.9409 + 2.0947 * VS - 0.8206 * VS ** 2.0 + 0.2683 * VS ** 3.0 - 0.0251 * VS ** 4.0'
#fld KEY VINF VSUP
#unt - - -
#fmt %5s %16f %16f
       -Z1        -0.333333        -0.001000
       -Z2        -0.666667        -0.333333
       -Z3        -1.000000        -0.666667
       -Z4        -1.333333        -1.000000
       -Z5        -1.666667        -1.333333
       -Z6        -2.000000        -1.666667
       -Z7        -2.333333        -2.000000
       -Z8        -2.666667        -2.333333
       -Z9        -3.000000        -2.666667
       VS0         1.000000         1.100000         
       VS1         1.100000         1.200000         
       VS2         1.200000         1.300000         
       VS3         1.300000         1.400000         
       VS4         1.400000         1.500000         
       VS5         1.500000         1.600000         
       VS6         1.600000         2.000000         
       VS7         2.000000         2.500000         
       VS8         2.500000         3.000000         
       VS9         3.000000         3.1
"""

    dispersion_file_string = """SURF96 R U T 0 5.000000 2.368112 0.1
SURF96 R U T 0 4.429334 2.289762 0.1
SURF96 R U T 0 3.923800 2.184182 0.1
SURF96 R U T 0 3.475964 2.028462 0.1
SURF96 R U T 0 3.079241 1.787782 0.1
SURF96 R U T 0 2.727797 1.457564 0.1
SURF96 R U T 0 2.416465 1.165715 0.1
SURF96 R U T 0 2.140666 1.051432 0.1
SURF96 R U T 0 1.896345 1.051072 0.1
SURF96 R U T 0 1.679909 1.085069 0.1
SURF96 R U T 0 1.488176 1.126465 0.1
SURF96 R U T 0 1.318325 1.168252 0.1
SURF96 R U T 0 1.167861 1.206051 0.1
SURF96 R U T 0 1.034569 1.234272 0.1
SURF96 R U T 0 0.916490 1.246419 0.1
SURF96 R U T 0 0.811888 1.235000 0.1
SURF96 R U T 0 0.719225 1.192478 0.1
SURF96 R U T 0 0.637137 1.115804 0.1
SURF96 R U T 0 0.564419 1.017582 0.1
SURF96 R U T 0 0.500000 0.929962 0.1
"""

    A = AsciiFile_fromstring(parameter_file_string)
    # p = Parameterizer_mZVSPRRH(A)
    p = Parameterizer_mZVSVPvsRHvp(A)
    d = makedatacoder(surf96filename=dispersion_file_string, which=Datacoder_log)

    t = Theory(parameterizer=p, datacoder=d)

    fd = t.frechet_derivatives(m=p.MMEAN, gm=None)

    import matplotlib.pyplot as plt
    plt.colorbar(plt.imshow(fd))
    plt.show()
