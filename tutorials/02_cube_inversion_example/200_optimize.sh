mkdir optimize
cd optimize || exit 1

## works
#HerrMet --optimize -init ../nodes.txt \
#  -data 1. 0.0 \
#  -prior  100. 1. \
#          200. 2. \
#          1 \
#          1.0 0.5 \
#  -upd 3 0 \
#  -save


# do not work
HerrMet --optimize -init ../nodes.txt \
  -data 1. 0.0 \
  -prior  30. 1. \
          40. 2. \
          1 \
          1.0 0.5 \
  -upd 3 0 \
  -save