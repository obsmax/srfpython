#!/usr/bin/env bash

#rm -f model000.10m.mod96
#m96 --split model000.mod96 -thck 0.01 -sfx 10m
#m96 --show  model000.50m.mod96 model000.mod96

sker17 -m96 model000.mod96 \
        -RU0 .2 2. 50 plog \
        -RU1 .2 2. 50 plog \
        -png
        # -norm  # use only for irregular depth sampling
