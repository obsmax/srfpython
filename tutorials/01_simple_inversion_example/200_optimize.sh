#!/usr/bin/env bash


# generate a template file
trash _HerrMet.optimize.param
HerrMet -verbose 0 --optimize -temp

# edit it and customize it ...
# in this tuto I edit here below :
cat << END > _HerrMet.optimize.param
# Parameter file for the optimization plugin
# ========== recall the grid parameters
nx           1
ny           1
dx           1.000000    # km
dy           1.000000    # km

# ========== path to the point-wise inversion data and results
# formattable strings including reference to iy and ix
datafiles    _HerrMet_data010/_HerrMet.target    
extractfiles _HerrMet_data010/_HerrMet.best_1000_0.00_1.p0.50.mod96

# ========== optimization parameters
ndecim       1  # downsamp the grid
newnz        100       # new vertical grid (needed to interpolate the extracted files)
newdz        0.03      # new vertical grid (needed to interpolate the extracted files)
Lh           1.000000   # horizontal smoothing distance km
Lv           0.100000   # vertical smoothing distance km
vsunc        0.100000   # km/s
damping      1.000000   # weight of the regularization

END

cat _HerrMet.optimize.param

