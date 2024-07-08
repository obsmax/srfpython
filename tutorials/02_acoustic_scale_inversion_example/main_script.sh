bash clear.sh

echo "# =========== control the picked dispersion curves"
echo "s96 --show picks.surf96 -freq"
s96 --show picks.surf96 -freq

echo "# =========== initiate an inversion dir using these picks as targert"
echo "HerrMet --target picks.surf96"
HerrMet --target picks.surf96

echo "# =========== Prepare 3 different inversions"
# cp -av _HerrMet_picks _HerrMet_4layers
# cp -av _HerrMet_picks _HerrMet_3layers
mv _HerrMet_picks _HerrMet_2layers

echo "# =========== Prepare 3 different parameterizations for comparison"
# Note : in most cases, you must edit the parameterization files before inversion
HerrMet --param 2 0.00004 -growing -t mZVSPRRH -op
HerrMet --send _HerrMet_2layers

# HerrMet --param 3 0.00004 -growing -t mZVSPRRH -op
# HerrMet --send _HerrMet_3layers

# HerrMet --param 4 0.00004 -growing -t mZVSPRRH -op
# HerrMet --send _HerrMet_4layers

echo "# =========== Check the target and prior boundaries"
echo "HerrMet --display"
HerrMet --display -si

echo "# =========== Run the inversion for all parameterizations"
echo "HerrMet --HerrMet --run"
HerrMet -w 12 --run -mode restart -nchain 4 -nkeep 1000

echo "# =========== Display results"
echo "HerrMet --HerrMet --display"
# -compact -si \
HerrMet --display \
    -plot best 500 0 1 \
    -pdf best 500 0 1 \
    -overdisp \
    -png
eog ./_HerrMet*/*png

