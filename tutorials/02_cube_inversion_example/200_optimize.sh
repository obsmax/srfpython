mkdir optimize
cd optimize || exit 1

# =============== using all nodes
## works
#HerrMet --optimize -init ../nodes.txt \
#  -data 1. 0.0 \
#  -prior  100. 1. \
#          200. 2. \
#          1 \
#          1.0 0.5 \
#  -upd 3 0 \
#  -save

# do not work if all nodes
HerrMet --optimize -init ../nodes.txt \
  -data 1. 0.0 \
  -prior  100. 1. \
          100. 2. \
          1 \
          1.0 0.5 \
  -upd 3 0 \
  -save

# still not work
#HerrMet --optimize -init ../nodes.txt \
#  -data 1. 0.0 \
#  -prior  100. 1. \
#          100. 2. \
#          0 \
#          1.0 0.5 \
#  -upd 3 0 \
#  -save

## do not work with some nodes
#HerrMet --optimize -init ../nodes.txt \
#  -data 1. 0.0 \
#  -prior  30. 1. \
#          30. 1. \
#          1 \
#          1.0 0.5 \
#  -upd 3 0 \
#  -save