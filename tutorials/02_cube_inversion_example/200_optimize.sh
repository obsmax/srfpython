mkdir optimize
cd optimize || exit 1

HerrMet --optimize -init ../nodes.txt
HerrMet --optimize -data -prior 100. 0. -fd

# HerrMet --optimize -prior 100. 0.
#exit 1
HerrMet --optimize -restart -upd 3  # -show

HerrMet --optimize -save