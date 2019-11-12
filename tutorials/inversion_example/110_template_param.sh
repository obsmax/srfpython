#!/usr/bin/env bash
HerrMet --param 9 3. -t mZVSPRRH -op \
    -dvp -.5 1.5 \
    -dvs -.5 1.5 \
    -drh -0. 1. \
    -dpr -1.0 0.

#--param      i f     generate a template parameter file to custom in .
#                     need the number of layers and bottom depth in km
#    -basedon s       build parametrization based on an existing mod96 file, require a filename, 
#                     if not specified, I take fixed values to build the parameter file
#    -t       s       parameterization type to use (['mZVSPRRH', 'mZVSVPRH', 'mZVSPRzRHvp', 'mZVSPRzRHz']), 
#                     default mZVSPRRH
#          mZVSPRRH = parameterize with 
#                     - depth interface (mZ1 = first interface, mZ2, ...), 
#                     - VS in each layer (VS0 = first layer, ...), 
#                     - VP/VS in each layer (PR0, PR1, ...), 
#                     - Density in each layer (RH0, RH1, ...)
#          mZVSVPRH = parameterize with  
#                     - depth interface (mZ1 = first interface, mZ2, ...), 
#                     - VS in each layer (VS0 = first layer, ...), 
#                     - VP in each layer (VP0, VP1, ...), 
#                     - Density in each layer (RH0, RH1, ...)
#       mZVSPRzRHvp = parameterize with  
#                     - depth interface (mZ1 = first interface, mZ2, ...), 
#                     - VS in each layer (VS0 = first layer, ...), 
#                     - use a fixed relation between VP/VS = f(z)
#                     - use a fixed relation between RH = f(VP)
#        mZVSPRzRHz = parameterize with  
#                     - depth interface (mZ1 = first interface, mZ2, ...), 
#                     - VS in each layer (VS0 = first layer, ...), 
#                     - use a fixed relation between VP/VS = f(z)
#                     - use a fixed relation between RH = f(z) 
#    -dvp     f f     add prior constraint on the vp offset between layers, 
#                     requires the extremal values, km/s
#    -dvs     f f     add prior constraint on the vs offset between layers, idem
#    -drh     f f     add prior constraint on the density offset between layers, idem, g/cm3
#    -dpr     f f     add prior constraint on the vp/vs offset between layers, idem, no unit
#    -growing         shortcut for -dvp 0. 5. -dvs 0. 5. -drh 0. 5. -dpr -5. 0.
#    -op              force overwriting ./_HerrMet.param if exists
#    
