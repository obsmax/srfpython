#!/usr/bin/env bash
#--run        s [s..] start inversion for the required rootnames, default _HerrMet_*
#    -mode    s       set the running mode, default skip
#                     restart  : overwrite the current run file(s) if any   
#                     append   : add new models to the exiting run file(s) 
#                     skip     : ignore rootnames with existsing run file(s)               
#    -nchain  i       number of chains to use, default 12
#    -nkeep   i       number of models to keep per chain, default 100
#    [use -w option before --run to control the maximum number of chains to run simultaneously]   

cd inversion || exit 1;

## quick test => very bad result
#HerrMet -w 8 --run  \
#    -mode restart \
#    -nchain 1 \
#    -nkeep 200
#exit 0;

HerrMet -w 8 --run  \
    -mode restart \
    -nchain 4 \
    -nkeep 1000
