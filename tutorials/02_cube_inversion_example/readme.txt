# erase previous test files
bash clear.sh

########## GENERATE A SYNTHETIC MODEL
# create a 3d depth model => .mod96
# the model is made of several "nodes" in the horizontal plane
python 000_create_3d_model.py

eog *png

# compute the dispersion curves at each node
bash 010_dispers.sh
bash 011_show.sh

########## INVERT

# create one directory per node for point-wise inversion
bash 100_set_target.sh
ls inversion

# define the parameters
# here I parameterize only vs, vp and rho are inferred from vs based on Brocher 2005
bash 111_param.sh

# see the parameterization
bash 112_show.sh

# copy the parameter file to the temporary directory of each node
# user may need to use a different parameterization for each node
bash 115_send.sh

# make sure the parameterization is ok
bash 117_show.sh

# run the inversion (~30mn on a 8 core machine) :
# custom 120_run.sh for the number of threads / markov chains to use
bash 120_run.sh

# see the chains status for all nodes, can be called while 120 is running
bash 130_manage.sh

# display the best models distributions => generate png files in each temp directory
bash 140_display.sh

# extract the inversion solution
bash 150_extract.sh

# compare inverted and expected vs velocity profiles
python 151_compare.py