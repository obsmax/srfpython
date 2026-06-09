#!/usr/bin/env bash

# compute dispersion curves from the command line
m96 --disp ./model000.mod96  \
  -RU0 .2 2. 20 plog \
  -RU1 .2 2. 20 plog  \
  -RC0 .2 2. 20 plog \
  -RC1 .2 2. 20 plog  \
  -save
  
# surf96 output has same name as mod96 file  by default=> rename it
mv model000.surf96 data010.surf96

