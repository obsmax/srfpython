#!/usr/bin/env bash

for n in 10 100 1000 10000
do
    HerrMet --display _HerrMet_data010 \
        -plot best $n 0 1 \
        -m96 model000.mod96 \
        -overdisp -png
        # -pdf  best $n 0 1 \

    mv _HerrMet_data010/_HerrMet.png  _HerrMet_data010/_HerrMet_best$n.png
done

eog _HerrMet_data010/_HerrMet_best*png

#HerrMet --display _HerrMet_data010 -plot best 1 0 1 -m96 model000.mod96
#HerrMet --display _HerrMet_data010 -plot best 100 0 1 -m96 model000.mod96 -overdisp
#HerrMet --display _HerrMet_data010 -plot best 1000 0 1 -m96 model000.mod96 -overdisp
#HerrMet --display _HerrMet_data010 -plot best 2000 0 1 -m96 model000.mod96 -overdisp -png
#HerrMet --display _HerrMet_data010 -plot best 0 -10 1 -m96 model000.mod96 -overdisp -png

#HerrMet --display _HerrMet_data010 -plot best 10000 0 1 -m96 model000.mod96  
#HerrMet --display _HerrMet_data010 -plot best 10000 0 1 -m96 model000.mod96  -png # -pdf best 10000 0 1 
#--display   s [s...] display param, target, and run outputs for the required rootnames, default _HerrMet_*
#                     (use "." to see the parameterzation template ./_HerrMet.param from option --param)
#    -plot  [s i f i] show the best models on the figure, arguments are :  
#                     first argument = selection mode, last or best
#                     second argument = highest model number to include (>=0, 0 means all)  
#                     third argument = lowest log likelyhood value to include (<=0.0, 0.0 means all)
#                     fourth argument = include only one model over "step" (>=1)
#                     default last 100 0.0 1             
#    -overdisp        recompute dispersion curves of the best models selected with higher resolution
#    -pdf   [s i f i] compute and show the statistics for the selected models, see -plot for arguments
#                     default last 0 0.0 1 
#                     use --extract to save pdf outputs
#    -png   [i]       save figure as pngfile instead of displaying it on screen, requires dpi, default 100 
#    -m96    s [s...] append depth model(s) to the plot from mod96 file(s)
#    -cmap            colormap, default viridis
#    -compact         display only vs and the dispersion curves, default False
#    -ftsz  i         set font size, default 10
#    -inline          do not pause (use in jupyter notebooks)
#    
