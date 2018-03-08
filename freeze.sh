`ps -ef | grep HerrMet | grep -v grep | awk 'BEGIN {printf("kill -SIGSTOP ")} {printf($2" ")} END {print " "}'`
