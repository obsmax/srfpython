#!/usr/bin/env bash


# generate a template file
trash _HerrMet.optimize.param
HerrMet -verbose 0 --optimize -temp

# edit it and customize it ...
# in this tuto I edit here below :
cat << END > _HerrMet.optimize.param
# Parameter file for the optimization plugin
# ========== Horizontal grid parameters
# assume that an inversion has run for each node in a cartesian grid defined by their nodes pairs (iy, ix)
# for a single inversion point use nx=ny=1, dx=dy=1.0
nx           5
ny           4
dx           2.4    # km
dy           2.2    # km

# ========== path to the point-wise inversion data and results
# use formattable strings including reference to iy and ix.

# Dispersion curve of each node, by default you should point to the _Herrmet*/_Herrmet.target
datafiles  ./inversion/_HerrMet_node_iy{iy:03d}_ix{ix:03d}/_HerrMet.target

# Which model to use as prior in each node (iy, ix)
# point either to the best fitting models (e.g. rank0) or to the median of n best or last models
# nb the start model is not the prior but the smooth version of it
priorfiles ./inversion/_HerrMet_node_iy{iy:03d}_ix{ix:03d}/_HerrMet.best_1000_0.00_1.p0.50.mod96

# Posterior uncertainty of McMC can be used as prior uncertainty of the optimization process
lowerfiles None
upperfiles None

# ========== optimization parameters
ndecim       1    # downsamp the grid, i.e. discard nodes of the horizontal grid
newnz        60   # new vertical grid (needed to interpolate the extracted files)
newdz        0.05 # new vertical grid (needed to interpolate the extracted files)
Lh           1.0  # horizontal smoothing distance km
Lv           0.5  # vertical smoothing distance km
vsunc        0.1  # velocity uncertainty km/s, or scaling factor to prior uncertainty if lower/upper files are provided
damping      1.0  # weight of the regularization
END

cat _HerrMet.optimize.param
