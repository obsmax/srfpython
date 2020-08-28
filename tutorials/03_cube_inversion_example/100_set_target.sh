#!/usr/bin/env bash

mkdir inversion
cd inversion || exit 1;

HerrMet --target ../data/node*.surf96 -lunc .1 -ot