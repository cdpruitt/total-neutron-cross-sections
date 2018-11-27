#!/bin/bash

expFile="/data2/analysis/total.root"
expGraph1Name="D_fromD2O"
expGraph2Name="H_fromH2O"

relGraph1Name="DtoH_exp"

litFile="/data2/analysis/literatureData.root"

litHGraphName="Hydrogen (n,tot)"
litDGraphName="Deuterium (n,tot)"

relGraph2Name="DtoH_lit"

outputFile="/data2/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expGraph1Name" "$expFile" "$expGraph2Name" "$outputFile" "$relGraph1Name"
../bin/relativeDiffCS "$litFile" "$litDGraphName" "$litFile" "$litHGraphName" "$outputFile" "$relGraph2Name"

../bin/multiplyCS "$outputFile" "$relGraph1Name" "100" "$relGraph1Name, percent"
../bin/multiplyCS "$outputFile" "$relGraph2Name" "100" "$relGraph2Name, percent"
