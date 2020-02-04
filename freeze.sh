#!/bin/bash

# a script to rapidely pause all running threads of HerrMet (using SIGSTOP)

`ps -ef | grep HerrMet | grep -v grep | awk 'BEGIN {printf("kill -SIGSTOP ")} {printf($2" ")} END {print " "}'`
