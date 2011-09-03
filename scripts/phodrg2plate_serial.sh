#!/bin/bash
# Takes the arguments:
# phodrg2plate <job file prefix> <url>

filename=$1$MPIEXEC_RANK
cat $filename | xargs -n1 phodrg2plate $2