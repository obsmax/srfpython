#!/bin/bash

# a script to restart all paused threads after freeze.sh

`ps -ef | grep HerrMet | grep -v grep | awk 'BEGIN {printf("kill -SIGCONT ")} {printf($2" ")} END {print " "}'`
