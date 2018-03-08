`ps -ef | grep HerrMet | grep -v grep | awk 'BEGIN {printf("kill -SIGCONT ")} {printf($2" ")} END {print " "}'`
