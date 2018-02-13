#!/bin/bash

expFile="/data2/analysis/corrected.root"
expGraph1Name="ONat_fromH2O"
expGraph2Name="O18"
expGraph3Name="ONat_fromD2O"

litFile="/data2/analysis/literatureData.root"
litGraph1Name="NatO(n,tot)"
litGraph2Name="18O(n,tot)"

relGraph1Name="ONat_fromH2O, expLit, corrected"
relGraph2Name="O18, expLit, corrected"
relGraph3Name="ONat_fromD2O, expLit, corrected"

outputFile="/data2/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expGraph1Name" "$litFile" "$litGraph1Name" "$outputFile" "$relGraph1Name"
../bin/relativeDiffCS "$expFile" "$expGraph2Name" "$litFile" "$litGraph2Name" "$outputFile" "$relGraph2Name"
../bin/relativeDiffCS "$expFile" "$expGraph3Name" "$litFile" "$litGraph1Name" "$outputFile" "$relGraph3Name"

# scale to percentage
../bin/multiplyCS "$outputFile" "$relGraph1Name" "100" "$relGraph1Name, percent"
../bin/multiplyCS "$outputFile" "$relGraph2Name" "100" "$relGraph2Name, percent"
../bin/multiplyCS "$outputFile" "$relGraph3Name" "100" "$relGraph3Name, percent"
