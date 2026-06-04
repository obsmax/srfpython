#!/usr/bin/env bash

# HerrMet --target -help
# --target     s [s..] A target designate a surf96 file containing a set of 
#                      dispersion curves to invert for one specific location.
#                      For each surf96 file, I create a directory in . for temporary files
#                      the target file in this directory is a copy of the input surf96 file. 
#                      It can be edited before inversion to remove unwanted points or resample dispersion curves... 
#                      it the target file is AAAA.surf96, the temporary directory will be named 
#                      ./_HerrMet_AAAA, this name is referred to as the "rootname" in the other plugins
#     -resamp  f f i s resample the dispersion curve in the target file, 
#                      requires fmin(Hz), fmax(Hz), nfreq, fscale 
#                      (flin=linear in freq domain, plin=linear in period or log=logarithmic scaling)
#     -lunc    f       set constant uncertainties in log domain (uncertainty = value x lunc)
#     -unc     f       set constant uncertainty in linear domain (uncertainty = unc)
#     -ot              force overwriting the targetfiles if exist
#     -h, -help        display the help message for this plugin


HerrMet --target data010.surf96 -lunc .1 -ot
#HerrMet --target data010.surf96 -unc .1 -ot
