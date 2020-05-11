mkdir optimize
cd optimize || exit 1

HerrMet --optimize -init ../nodes.txt

HerrMet --optimize \
  -data 1. 0.0

HerrMet --optimize \
  -prior  100. 1. \
          200. 2. \
          1 \
          1.0 0.5

# run two iterations without updating the fd
HerrMet --optimize -restart -upd 2 0
# update the fd manually
HerrMet --optimize -fd
# run one more iteration with the new fd
HerrMet --optimize -upd 1 0

HerrMet --optimize -save