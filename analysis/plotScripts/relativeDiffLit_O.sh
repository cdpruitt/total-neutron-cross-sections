#!/bin/bash

expFile="/data2/analysis/total.root"
expGraph1Name="ONat_fromH2O"
expGraph2Name="O18"
expGraph3Name="ONat_fromD2O"

expGraphSE1Name="ONat_fromH2OSysErrors"
expGraphSE2Name="O18SysErrors"
expGraphSE3Name="ONat_fromD2OSysErrors"

litFile="/data2/analysis/literatureData.root"
litGraph1Name="NatO(n,tot)"
litGraph2Name="18O(n,tot)"

relGraph1Name="ONat_fromH2O, expLit"
relGraph2Name="O18, expLit"
relGraph3Name="ONat_fromD2O, expLit"

relGraphSE1Name="ONat_fromH2O, expLit, SysErrors"
relGraphSE2Name="O18, expLit, SysErrors"
relGraphSE3Name="ONat_fromD2O, expLit, SysErrors"

outputFile="/data2/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expGraph1Name" "$litFile" "$litGraph1Name" "$outputFile" "$relGraph1Name"
../bin/relativeDiffCS "$expFile" "$expGraph2Name" "$litFile" "$litGraph2Name" "$outputFile" "$relGraph2Name"
../bin/relativeDiffCS "$expFile" "$expGraph3Name" "$litFile" "$litGraph1Name" "$outputFile" "$relGraph3Name"

# scale to percentage
../bin/multiplyCS "$outputFile" "$relGraph1Name" "100" "$relGraph1Name, percent"
../bin/multiplyCS "$outputFile" "$relGraph2Name" "100" "$relGraph2Name, percent"
../bin/multiplyCS "$outputFile" "$relGraph3Name" "100" "$relGraph3Name, percent"

# calculate for systematic errors graphs
../bin/relativeDiffCS "$expFile" "$expGraphSE1Name" "$litFile" "$litGraph1Name" "$outputFile" "$relGraphSE1Name"
../bin/relativeDiffCS "$expFile" "$expGraphSE2Name" "$litFile" "$litGraph2Name" "$outputFile" "$relGraphSE2Name"
../bin/relativeDiffCS "$expFile" "$expGraphSE3Name" "$litFile" "$litGraph1Name" "$outputFile" "$relGraphSE3Name"

# scale to percentage
../bin/multiplyCS "$outputFile" "$relGraphSE1Name" "100" "$relGraphSE1Name, percent"
../bin/multiplyCS "$outputFile" "$relGraphSE2Name" "100" "$relGraphSE2Name, percent"
../bin/multiplyCS "$outputFile" "$relGraphSE3Name" "100" "$relGraphSE3Name, percent"
