#!/bin/bash

expFile="/data2/analysis/total.root"

inputGraph1Name="H2O"
inputGraph2Name="H2O18"

outputGraph1Name="ONat"
outputGraph2Name="O18"

litFile="/data2/analysis/literatureData.root"

litHGraphName="Hydrogen (n,tot)"

../bin/subtractCS "$expFile" "$inputGraph1Name" "$litFile" "$litHGraphName" 2 1 "$outputGraph1Name"
../bin/subtractCS "$expFile" "$inputGraph2Name" "$litFile" "$litHGraphName" 2 1 "$outputGraph2Name"

# shift O18 cross section up by 1 barn, for readability on plots

shiftFileName="/data2/analysis/shifted.root"
shiftGraphName="$outputGraph2Name+1Barn"

../bin/shiftCS "$expFile" "$outputGraph2Name" "1" "$shiftFileName" "$shiftGraphName"
