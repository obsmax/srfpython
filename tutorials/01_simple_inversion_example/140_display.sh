#!/usr/bin/env bash

# HerrMet --display -help
# --display   s [s...] display param, target, and run outputs for the required rootnames, default ./_HerrMet_*
#                      (use "." to see the parameterzation template _HerrMet.param from option --param)
#     -plot  [s i f i] show the best models on the figure, arguments are :  
#                      first argument = selection mode, last or best
#                      second argument = highest model number to include (>=0, 0 means all)  
#                      third argument = lowest log likelihood value to include (<=0.0, 0.0 means all)
#                      fourth argument = include only one model over "step" (>=1)
#                      default last 100 0.0 1             
#     -overdisp        recompute dispersion curves of the best models selected with higher resolution
#     -ddc             solver accuracy to be used if overdisp is used, it is adviced to use the same values as the run, default 0.005
#     -hstep           solver step to convert phase to group if overdisp is used, it is adviced to use the same values as the run, default 0.005    
#     -overmode i      use together with overdisp to add more overtones [experimental], 
#                      provide max mode number
#     -pdf   [s i f i] compute and show the statistics for the selected models, see -plot for arguments
#                      default last 0 0.0 1 
#                      use --extract to save pdf outputs
#     -png/-svg [i]    save figure instead of displaying it on screen, requires dpi, default 100 
#     -m96    s [s...] append depth model(s) to the plot from mod96 file(s)
#     -cmap            colormap, default viridis
#     -compact         display only vs and the dispersion curves, default False
#     -si              display in units of the international system
#     -ftsz  i         set font size, default 10
#     -inline          do not pause (use in jupyter notebooks)
#     -h, -help        display the help message for this plugin  
    
HerrMet --display _HerrMet_data010 \
    -plot best 1000 0 1 \
    -pdf  best 1000 0 1 \
    -ddc 0.002 \
    -compact -si \
    -m96 model000.mod96 \
    -overdisp  


