echo -e "\n# erase previous test files "
echo "    bash clear.sh ..."
bash clear.sh


echo -e "\n########## GENERATE A SYNTHETIC MODEL"
echo "    # create a depth model => .mod96"
echo "        python3 000_create_model.py ..."
python3 000_create_model.py

echo -e "\n    # show it"
echo "        bash 001_show.sh ..."
bash 001_show.sh

echo -e "\n    # compute the sensitivity kernels (optional)"
echo "        bash 002_sensitivity.sh ..."
bash 002_sensitivity.sh

echo -e "\n    # compute the dispersion curves"
echo "        bash 010_dispers.sh ..."
bash 010_dispers.sh

echo -e "\n    # show it"
echo "        bash 011_show.sh ..."
bash 011_show.sh


echo "########## INVERSE THE SYNTHETIC DATA"
echo -e "\n    # set the forward dispersion curves as the inversion data target"
echo "    # => create a subdirectory for each target disp. curve  (_HerrMet_*)"
echo "        bash 100_set_target.sh ..."
bash 100_set_target.sh

echo -e "\n    # show the current state of the inversion in the sub directory : only data for now"
echo "        bash 101_show.sh ..."
bash 101_show.sh

echo -e "\n    # build a template parameter file using HerrMet in ."
echo "        bash 110_template_param.sh ..."
bash 110_template_param.sh

echo -e "\n    # simulate the user editing the parameter file, here I simply overwrite the template using cat"
echo "        bash 111_param.sh ..."
bash 111_param.sh

echo -e "\n    # show the content of . : only parameterization since the targets are in the sub directories (_HerrMet_*)"
echo "        bash 112_show.sh ..."
bash 112_show.sh

echo -e "\n    # send the custom parameter file from . to all the subdirectories"
echo "    # here there is only one"
echo "        bash 115_send.sh ..."
bash 115_send.sh

echo -e "\n    # show the current state of the inversion in the sub directory : data and parameterization"
echo "        bash 117_show.sh ..."
bash 117_show.sh

echo -e "\n    # run the inversion in parallel in the subdirectories"
echo "        bash 120_run.sh ..."
bash 120_run.sh

echo -e "\n    # show the statistics of the inversion"
echo "        bash 130_manage.sh ..."
bash 130_manage.sh

echo -e "\n    # create figures"
echo "        bash 140_display.sh ..."
bash 140_display.sh

echo -e "\n    # extract the inversion solution for later use..."
echo "        bash 150_extract.sh ..."
bash 150_extract.sh

