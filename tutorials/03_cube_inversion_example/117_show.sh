#!/usr/bin/env bash

cd inversion || exit 1;

HerrMet --display  _HerrMet_node_iy007_ix001 -m96 ../model/nodes/node_iy007_ix001.mod96 -png

eog _HerrMet_node*/_HerrMet.png