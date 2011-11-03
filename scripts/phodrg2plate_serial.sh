#!/bin/bash
# Takes the arguments:
# phodrg2plate <job file prefix> <url>

filename=$1`echo $MPIEXEC_RANK | awk '{printf("%.2d\n", $1);}'`
cat $filename | xargs -n1 phodrg2plate $2