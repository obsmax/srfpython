mkdir optimize
cd optimize || exit 1

HerrMet --optimize -init ../nodes.txt
## WORKS
#HerrMet --optimize -data \
#  -prior  100. 1. 200. 2. 1 1.0 0.1 \
#  -fd

## WORKS
#HerrMet --optimize -data \
#  -prior  100. 2. 200. 4. 1 1.0 0.1 \
#  -fd

HerrMet --optimize -data \
  -prior  100. 1. 200. 2. 1 1.0 0.0 \
  -fd

# HerrMet --optimize -prior 100. 0.
#exit 1
HerrMet --optimize -restart -upd 1  -show

HerrMet --optimize -save