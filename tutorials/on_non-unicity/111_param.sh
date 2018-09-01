cat << END > _HerrMet.param
#met NLAYER = 5
#met TYPE = 'mZVSPRRH'
#fld KEY VINF VSUP
#unt - - -
#fmt %5s %16f %16f
       -Z1        -0.5             -0.05  
       -Z2        -1.5             -0.25
       -Z3        -2.5             -1.0
       -Z4        -3.0             -1.2
       VS0         0.100000         3.100000
       VS1         0.200000         3.200000
       VS2         0.300000         3.300000
       VS3         0.400000         3.400000
       VS4         0.50000          3.500000
       PR0         1.90             2.1
       PR1         1.80             2.0
       PR2         1.75             1.9
       PR3         1.71             1.85
       PR4         1.70             1.8
       RH0         2.10             2.5
       RH1         2.20             2.6
       RH2         2.30             2.7
       RH3         2.40             2.7
       RH4         2.50             2.7
END
