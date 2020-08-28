#!/usr/bin/env bash
cd inversion || exit 1;

cat << END > _HerrMet.param
#met  NLAYER = 5
#met TYPE = 'mZVSVPvsRHvp'
#met PRIORTYPE = 'DVPDVSDRH'
#met DVSMIN = -0.5
#met DVSMAX = +1.5
#met DVPMIN = -0.5
#met DVPMAX = +1.5
#met DPRMIN = -1.0
#met DPRMAX = +0.0
#met DRHMIN = +0.0
#met DRHMAX = +1.0

#met VPvs = 'lambda VS: 0.9409 + 2.0947 * VS - 0.8206 * VS ** 2.0 + 0.2683 * VS ** 3.0 - 0.0251 * VS ** 4.0'
#met RHvp = 'lambda VP: 1.6612 * VP - 0.4721 * VP ** 2.0 + 0.0671 * VP ** 3.0 - 0.0043 * VP ** 4.0 + 0.000106 * VP ** 5.0'

#fld KEY     VINF          VSUP
#unt []      []            []
#fmt %s      %f            %f
     -Z1     -0.5          -0.01
     -Z2     -1.0          -0.25
     -Z3     -2.0          -0.75
     -Z4     -3.0          -1.5

     VS0     0.2           2.0
     VS1     0.5           2.5
     VS2     1.0           3.0
     VS3     1.5           3.25
     VS4     2.0           3.9


END
