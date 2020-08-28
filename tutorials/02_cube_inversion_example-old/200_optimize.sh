mkdir optimize
cd optimize || exit 1

# =============== using all nodes
HerrMet --optimize -init ../nodes.txt \
  -data 1. 0.0 \
  -prior  100. 2.0 \
          1000. 2000.0 \
          0 \
          0.0 0.8 \
  -upd 5 0 \
  -save
