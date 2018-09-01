rm -f model000.50m.mod96
m96 --split model000.mod96 -thck 0.05 -sfx 50m

#m96 --show  model000.50m.mod96 model000.mod96

sker17 -m96 model000.50m.mod96 \
        -RU0 .2 1. 50 plog \
        -RC0 .2 1. 50 plog \
        -norm
