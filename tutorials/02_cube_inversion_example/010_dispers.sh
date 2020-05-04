#!/usr/bin/env bash

for f in models/node???.mod96;
do
  m96 --disp $f  \
        -RU0 .2 2. 20 plog \
        -RU1 .2 2. 20 plog  -save
done

mkdir --parents data
rm -f data/node???.surf96

mv ./node???.surf96 data/
