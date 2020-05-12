mkdir optimize
cd optimize || exit 1

HerrMet --optimize -init ../nodes.txt \
  -data 1. 0.0 \
  -prior  100. 1. \
          200. 2. \
          1 \
          1.0 0.5 \
  -upd 3 0 \
  -save