#!/usr/bin/env bash

cd inversion || exit 1;

HerrMet --display  _HerrMet_node_iy002_ix000 -m96 ../model/nodes/node_iy002_ix000.mod96 -png

eog _HerrMet_node*/_HerrMet.png