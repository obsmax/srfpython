#!/usr/bin/env bash


# generate a template file
trash _HerrMet.optimize.param
HerrMet -verbose 0 --optimize -temp

# edit it and customize it ...
# in this tuto I edit here below :
cat << END > _HerrMet.optimize.param
# ========== recall the grid parameters
nx           15
ny           10
dx           2.4    # km
dy           2.2    # km

# ========== path to the point-wise inversion data and results
datafiles    ./inversion/_HerrMet_node_iy{iy:03d}_ix{ix:03d}/_HerrMet.target
extractfiles ./inversion/_HerrMet_node_iy{iy:03d}_ix{ix:03d}/_HerrMet.best_1000_0_1.p0.50.mod96

# ========== optimization parameters
ndecim       1    # downsamp the grid
newnz        10   # new vertical grid (needed to interpolate the extracted files)
newdz        0.33 # new vertical grid (needed to interpolate the extracted files)
Lh           1.0  # horizontal smoothing distance km
Lv           0.5  # vertical smoothing distance km
vsunc        0.1  # km/s
damping      1.0  # weight of the regularization
END

cat _HerrMet.optimize.param
