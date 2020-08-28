#!/usr/bin/env bash

cd inversion || exit 1;
for node in 000 001 002 003 004 005 006 007 008 009;
do
  HerrMet --display _HerrMet_node$node -m96  ../models/node$node.mod96 -png
done
eog _HerrMet_node00?/_HerrMet.png