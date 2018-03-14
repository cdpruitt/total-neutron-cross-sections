#!/bin/bash

expFile="/data2/analysis/total.root"
expGraph1Name="CNat"
expGraph2Name="shortCNat"

relGraph1Name="longCShortC"

expGraphSE1Name="CNatSysErrors"
expGraphSE2Name="shortCNatSysErrors"

relGraphSE1Name="longCShortCSysErrors"

outputFile="/data2/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expGraph1Name" "$expFile" "$expGraph2Name" "$outputFile" "$relGraph1Name"

../bin/multiplyCS "$outputFile" "$relGraph1Name" "100" "$relGraph1Name, percent"

# sys errors
../bin/relativeDiffCS "$expFile" "$expGraphSE1Name" "$expFile" "$expGraphSE2Name" "$outputFile" "$relGraphSE1Name"

../bin/multiplyCS "$outputFile" "$relGraphSE1Name" "100" "$relGraphSE1Name, percent"
