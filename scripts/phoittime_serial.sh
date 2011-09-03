#!/bin/bash
# Takes the arguments:
# phoitalbedo_serial.sh <number of jobs> <level> <url>
phoittime -l $2 -j $MPIEXEC_RANK -n $1 $3