# erase previous test files
bash clear.sh


########## GENERATE A SYNTHETIC MODEL TO INVERSE
# create a depth model => .mod96
python 000_create_model.py

# show it
bash 001_show.sh

# compute the sensitivity kernels (optional)
bash 002_sensitivity.sh

# compute the dispersion curves
bash 010_dispers.sh

# show it
bash 011_show.sh


########## INVERSE THE SYNTHETIC DATA
# set the forward dispersion curves as the inversion data target
# => create a subdirectory for each target disp. curve  (_HerrMet_*)
bash 100_set_target.sh

# show the current state of the inversion in the sub directory : only data for now
bash 101_show.sh

# build a template parameter file using HerrMet in .
bash 110_template_param.sh

# simulate the user editing the parameter file, here I simply overwrite the template using cat
bash 111_param.sh

# show the content of . : only parameterization since the targets are in the sub directories (_HerrMet_*)
bash 112_show.sh

# send the custom parameter file from . to all the subdirectories
# here there is only one
bash 115_send.sh

# show the current state of the inversion in the sub directory : data and parameterization
bash 117_show.sh

# run the inversion in parallel in the subdirectories
bash 120_run.sh

# show the statistics of the inversion
bash 130_manage.sh

# create figures
bash 140_display.sh

# extract the inversion solution for later use...
bash 150_extract.sh

