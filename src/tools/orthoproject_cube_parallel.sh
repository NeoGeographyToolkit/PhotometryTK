#!/bin/bash

ORTHOPROJECT=$1
tile=$2
cub=$3
adj=$4
mpp=$5
outDir=$6

tilePrefix=$(echo $tile | perl -pi -e "s#^.*\/##g" | perl -pi -e "s#\..*?\$##g")

outputImg="$outDir/$tilePrefix.tif"

$ORTHOPROJECT --mark-no-processed-data --mpp $mpp $tile $cub $adj $outputImg


