#!/usr/bin/env bash
cd inversion || exit 1;

# remove previous extracted files if any
find inversion -name "_HerrMet.best_*.p0.??.*96" -delete

HerrMet --extract  \
    -pdf best 1000 0 1
