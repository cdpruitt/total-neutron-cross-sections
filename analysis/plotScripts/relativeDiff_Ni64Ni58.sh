#!/bin/bash

expFile="/data1/analysis/total.root"
expGraph1Name="Ni64PurityCorrected"
expGraph2Name="Ni58"

relGraph1Name="Ni64Ni58"

expGraphSE1Name="Ni64PurityCorrectedSysErrors"
expGraphSE2Name="Ni58SysErrors"

relGraphSE1Name="Ni64Ni58SysErrors"

outputFile="/data1/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expGraph1Name" "$expFile" "$expGraph2Name" "$outputFile" "$relGraph1Name"

../bin/multiplyCS "$outputFile" "$relGraph1Name" "100" "$relGraph1Name, percent"

# systematic errors graphs
../bin/relativeDiffCS "$expFile" "$expGraphSE1Name" "$expFile" "$expGraphSE2Name" "$outputFile" "$relGraphSE1Name"

../bin/multiplyCS "$outputFile" "$relGraphSE1Name" "100" "$relGraphSE1Name, percent"
