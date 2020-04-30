#!/usr/bin/env bash
cd inversion || exit 1;

# I use a loop to compare each output with the expected model
for rootname in _HerrMet_node???
do
  node=`awk 'BEGIN {split("'$rootname'",a,"_");print a[3]}'`
  echo $node
  HerrMet --display  $rootname \
      -plot best 1000 0 1 \
      -pdf best 1000 0 1 \
      -m96 ../models/$node.mod96 \
      -png \
      || exit 1 ;
#      -overdisp \
done
eog _HerrMet_*/_HerrMet.png

exit 0;

# If you don't want a m96 file to appear on each plot, the following command is ok
HerrMet --display   \
      -plot best 1000 0 1 \
      -pdf best 1000 0 1 \
      -png
