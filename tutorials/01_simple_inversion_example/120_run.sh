#!/usr/bin/env bash

# 

# HerrMet --run -help
# --run        s [s..] start inversion for the required rootnames, default ./_HerrMet_*
#     -mode    s       set the running mode, default skip
#                      restart  : overwrite the current run file(s) if any   
#                      append   : add new models to the exiting run file(s) 
#                      skip     : ignore rootnames with existsing run file(s)               
#     -nchain  i       number of chains to use, default 12
#     -nkeep   i       number of models to keep per chain, default 100
#     -surflim f  f    Controls the uper and lower surface wave velocity
#                      default 0.025, 4.8km/s
#     -notarget        replaces the data term of the cost functions
#                      by a n-dimensional uniform PDF in the dispersion space
#                      between the two limits provided by -surflim
#                      This allows to visualize the dataspace mapped 
#                      according to the parameterization and the prior conditions                         
#     -ddc             solver accuracy, default 0.005
#     -timeout f       time limit for the Markov Chains in seconds, default None  
#     -hstep           solver step to convert phase to group, default 0.005
#     -h, -help        display the help message for this plugin
#     [use -w option before --run to control the maximum number of chains to run simultaneously]
  
    
HerrMet -w 8 --run _HerrMet_data010 \
    -mode restart \
    -ddc 0.002 \
    -nchain 8 \
    -nkeep 2000 \
    -timeout 120

# The output prompt gives you a status of the run progression
# the following line means    
# data010    chain    6 kept   34/  400 elapsed    5/   60s fail  114 stay   21 MP  0.13 AK 0.09 IK 0.00 AS 84.56/s IS 63.82/s BST -2349.285358

#  root  folder ./_HerrMet_data010
#  chain number 6
#  has   kept 34 models over 400 tested
#  has   run for 5 seconds over the 60s limit
#  has   failed solving the forward problem 114 times (inconsistent models tested)
#  stay  is the number of subsequent times a Markov Chain fails at moving, which weights the current model
#  MP    master proposal : is a coefficient to adjust the exploration area, 
#        this coeff adjusts automatically to maintain an average keep ratio at ~25%
#  AK    Average Keep ratio, the ratio of models kept over trashed
#  IK    Instantaneous Keep ratio, same as AK but for recent search
#  AS    Average speed, the number of model tested per second on average
#  IS    Instantaneous speed, the current number of models tested per second for this chain
#  BST   the log likelihood of the best model found for this chain






# 

