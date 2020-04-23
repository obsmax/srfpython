#!/usr/bin/env bash

rm -f ./_HerrMet_data010/*rank*mod96
rm -f ./_HerrMet_data010/*best_*mod96

# extract the 10 best models and 16% 50% 84% percentiles of the 10^4 best models
HerrMet --extract _HerrMet_data010 \
    -top 10 0 1 \
    -pdf best 10000 0 1

m96 --show \
    _HerrMet_data010/*rank*.mod96 \
    _HerrMet_data010/_HerrMet.best_10000_0_1.p*.mod96
