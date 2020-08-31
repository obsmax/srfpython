#!/usr/bin/env bash


# run the interative inversion
HerrMet --optimize -inv

# show the results
HerrMet --optimize -show y 10
