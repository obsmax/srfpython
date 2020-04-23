HerrMet --manage _HerrMet_data010 -stats -plot -100.

#--manage     s [s..] manage run results for given rootnames, default _HerrMet_*
#     -stats          prints detailed stats for each chain of each runfile 
#     -plot   [f]     display the convergence for every chain and every rootname, specify the lower bound
#     -inline         do not pause (jupyter)
#     -delbad f       delete bad models, log likelihood below a given threshold, no default
#     -delchains i [i...] delete one or more chains using their chainid
