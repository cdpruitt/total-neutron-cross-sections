#!/bin/bash

expFile="/data2/analysis/total.root"
expGraph1Name="Sn124"
expGraph2Name="Sn112"

relGraph1Name="Sn124Sn112"

expGraphSE1Name="Sn124SysErrors"
expGraphSE2Name="Sn112SysErrors"

relGraphSE1Name="Sn124Sn112SysErrors"

outputFile="/data2/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expGraph1Name" "$expFile" "$expGraph2Name" "$outputFile" "$relGraph1Name"

../bin/multiplyCS "$outputFile" "$relGraph1Name" "100" "$relGraph1Name, percent"

# systematic errors graphs
../bin/relativeDiffCS "$expFile" "$expGraphSE1Name" "$expFile" "$expGraphSE2Name" "$outputFile" "$relGraphSE1Name"

../bin/multiplyCS "$outputFile" "$relGraphSE1Name" "100" "$relGraphSE1Name, percent"
