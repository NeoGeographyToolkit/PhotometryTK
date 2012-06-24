#!/bin/bash
# __BEGIN_LICENSE__
#  Copyright (c) 2009-2012, United States Government as represented by the
#  Administrator of the National Aeronautics and Space Administration. All
#  rights reserved.
#
#  The NGT platform is licensed under the Apache License, Version 2.0 (the
#  "License"); you may not use this file except in compliance with the
#  License. You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
# __END_LICENSE__

ORTHOPROJECT=$1
tile=$2
cub=$3
adj=$4
mpp=$5
outDir=$6

tilePrefix=$(echo $tile | perl -pi -e "s#^.*\/##g" | perl -pi -e "s#\..*?\$##g")

outputImg="$outDir/$tilePrefix.tif"

$ORTHOPROJECT --mark-no-processed-data --mpp $mpp $tile $cub $adj $outputImg


