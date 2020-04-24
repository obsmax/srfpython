#!/usr/bin/env bash

cat << END > _HerrMet.param
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
END



#cat << END > _HerrMet.param
##met NLAYER = 5
##met TYPE = 'mZVSPRRH'
##fld KEY VINF VSUP
##unt - - -
##fmt %5s %16f %16f
#       -Z1        -0.2             -0.05
#       -Z2        -1.0             -0.10
#       -Z3        -2.5             -0.5
#       -Z4        -3.0             -1.2
#       VS0         0.100000         3.100000
#       VS1         0.200000         3.200000
#       VS2         0.300000         3.300000
#       VS3         0.400000         3.400000
#       VS4         0.50000          3.500000
#       PR0         1.90             2.1
#       PR1         1.80             2.0
#       PR2         1.75             1.9
#       PR3         1.71             1.85
#       PR4         1.70             1.8
#       RH0         1.7              2.5
#       RH1         2.0              2.6
#       RH2         2.20             2.7
#       RH3         2.40             2.7
#       RH4         2.50             2.7
#END
