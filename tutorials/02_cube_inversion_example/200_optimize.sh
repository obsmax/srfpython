mkdir optimize
cd optimize || exit 1

HerrMet --optimize -init ../nodes.txt

HerrMet --optimize \
  -data \
  -prior  100. 1. \
          200. 2. \
          1 \
          1.0 0.5

HerrMet --optimize -restart -upd 3 #  -show

HerrMet --optimize -save