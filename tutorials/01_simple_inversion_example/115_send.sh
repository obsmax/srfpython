#!/bin/bash

# copy the current parameterization file to 
# the desired inversion runs (can be all or some specific inversions)

#HerrMet --send -help
#--send       s [s..] send the custom parameterization file _HerrMet.param to the specified rootnames, 
#                     default ./_HerrMet_*
#    -op              force overwriting [rootname]/_HerrMet.param if exists
#    -h, -help        display the help message for this plugin


HerrMet --send -op
