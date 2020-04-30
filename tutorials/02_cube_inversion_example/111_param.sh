#!/usr/bin/env bash
cd inversion || exit 1;

cat << END > _HerrMet.param
#met  NLAYER = 9
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
     -Z1     -0.25         -0.01
     -Z2     -0.5          -0.1
     -Z3     -0.8          -0.2
     -Z4     -1.0          -0.3
     -Z5     -1.5          -0.5
     -Z6     -2.0          -0.8
     -Z7     -2.5          -1.0
     -Z8     -3.0          -1.5

     VS0     0.2           2.0
     VS1     0.5           2.1
     VS2     0.5           2.2
     VS3     0.6           2.3
     VS4     0.7           2.5
     VS5     1.0           2.8
     VS6     1.5           3.0
     VS7     2.0           3.2
     VS8     2.2           3.8

END
