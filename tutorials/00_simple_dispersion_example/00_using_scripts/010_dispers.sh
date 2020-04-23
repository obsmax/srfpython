#!/usr/bin/env bash

m96 --disp ./model000.mod96  -RU0 .2 2. 20 plog -RU1 .2 2. 20 plog  -save
mv model000.surf96 data010.surf96

