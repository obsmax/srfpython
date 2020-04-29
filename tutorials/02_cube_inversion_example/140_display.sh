#!/usr/bin/env bash
cd inversion || exit 1;

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
