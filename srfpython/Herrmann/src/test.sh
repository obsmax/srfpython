#!/bin/bash

#test the modified version of fortran programs max_srfpre96 and max_srfdis96


rm -f in.txt out1.txt out2.txt

# write input arguments to pass to srf_pre96 through stdin
cat << END > in1.txt
2 2 2 2
SURF96 L C T 0 7.692308 3.152852 0.1
SURF96 L C T 0 4.479947 2.655829 0.1
SURF96 L C T 0 2.609090 1.524800 0.1
SURF96 L C T 0 1.519516 1.129654 0.1
SURF96 L C T 0 0.884956 0.983793 0.1
SURF96 L C T 1 1.605682 3.069108 0.1
SURF96 L C T 1 0.900901 1.713844 0.1
SURF96 L U T 0 12.500000 3.148224 0.1
SURF96 L U T 0 6.681531 2.618293 0.1
SURF96 L U T 0 3.571429 0.983935 0.1
SURF96 L U T 0 1.909009 0.831358 0.1
SURF96 L U T 0 1.020408 0.832491 0.1
SURF96 L U T 1 1.839556 2.590743 0.1
SURF96 L U T 1 1.010101 0.977466 0.1
SURF96 R C T 0 6.666667 2.800877 0.1
SURF96 R C T 0 4.006426 2.617169 0.1
SURF96 R C T 0 2.407717 2.268805 0.1
SURF96 R C T 0 1.446951 1.284619 0.1
SURF96 R C T 0 0.869565 0.968552 0.1
SURF96 R C T 1 2.242243 2.193912 0.1
SURF96 R C T 1 1.384358 1.237185 0.1
SURF96 R C T 1 0.854701 0.962565 0.1
SURF96 R U T 0 10.000000 2.727757 0.1
SURF96 R U T 0 5.623413 2.463522 0.1
SURF96 R U T 0 3.162278 1.964674 0.1
SURF96 R U T 0 1.778279 0.642590 0.1
SURF96 R U T 0 1.000000 0.714381 0.1
SURF96 R U T 1 2.041241 0.947028 0.1
SURF96 R U T 1 1.304237 0.696091 0.1
SURF96 R U T 1 0.833333 0.708612 0.1
END

# call max_srfpre96, save outputs to be compared with expected ones
../bin/max_srfpre96 < in1.txt > out1.txt

# write input arguments to be passed to max_srfdis96, save outputs ...
cat << END > in2.txt
  8
   0.250   0.200   0.200   0.200   0.200   0.480   0.270
  1.850  2.360  2.630  3.150  3.710  4.540  5.480  5.800
  0.860  1.100  1.240  1.470  1.730  2.130  3.130  3.310
 2.470 2.470 2.470 2.470 2.470 2.580 2.580 2.630
  0.005 0.005
END
cat out1.txt >> in2.txt

../bin/max_srfdis96 < in2.txt > out2.txt

# compare output with expectations
if cmp --silent out1.txt expected1.txt; then
    echo "test1 ok"
else
  echo "error : the output from max_srfpre96 (out1.txt) differs from expected1.txt" \
  exit 1;
fi

if cmp --silent out2.txt expected2.txt; then
    echo "test2 ok"
else
  echo "error : the output from max_srfdis96 (out2.txt) differs from expected2.txt"; \
  exit 1;
fi

rm -f in.txt out1.txt out2.txt
exit 0;
