#!/bin/bash

expFile="/data2/analysis/corrected.root"

inputGraph1Name="H2O"
inputGraph2Name="H2O18"
inputGraph3Name="D2O"

outputGraph1Name="ONat_fromH2O"
outputGraph2Name="O18"
outputGraph3Name="ONat_fromD2O"

litFile="/data2/analysis/literatureData.root"

litHGraphName="Hydrogen (n,tot)"
litDGraphName="Deuterium (n,tot)"

../bin/subtractCS "$expFile" "$inputGraph1Name" "$litFile" "$litHGraphName" 2 1 "$outputGraph1Name"
../bin/subtractCS "$expFile" "$inputGraph2Name" "$litFile" "$litHGraphName" 2 1 "$outputGraph2Name"
../bin/subtractCS "$expFile" "$inputGraph3Name" "$litFile" "$litDGraphName" 2 1 "$outputGraph3Name"

# shift O18 cross section up by 1 barn, for readability on plots

shiftFileName="/data2/analysis/shifted_Corrected.root"
shiftGraphName="$outputGraph2Name+1Barn"

../bin/shiftCS "$expFile" "$outputGraph2Name" "1" "$shiftFileName" "$shiftGraphName"
