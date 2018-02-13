#!/bin/bash

expFile="/data2/analysis/total.root"
expGraph1Name="O18"
expGraph2Name="ONat_fromH2O"

relGraph1Name="O18O16"

outputFile="/data2/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expGraph1Name" "$expFile" "$expGraph2Name" "$outputFile" "$relGraph1Name"

../bin/multiplyCS "$outputFile" "$relGraph1Name" "100" "$relGraph1Name, percent"
