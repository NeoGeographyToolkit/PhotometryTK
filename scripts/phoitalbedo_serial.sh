#!/bin/bash
# Takes the arguments:
# phoitalbedo_serial.sh <number of jobs> <level> <url>
phoitalbedo -l $2 -j $MPIEXEC_RANK -n $[$1*3] $3
phoitalbedo -l $2 -j $[$MPIEXEC_RANK+$1] -n $[$1*3] $3
phoitalbedo -l $2 -j $[$MPIEXEC_RANK+$1+$1] -n $[$1*3] $3