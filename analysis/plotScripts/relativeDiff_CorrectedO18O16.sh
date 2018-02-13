#!/bin/bash

expFile="/data2/analysis/corrected.root"
expGraph1Name="O18"
expGraph2Name="ONat_fromH2O"

relGraph1Name="O18O16_corrected"

outputFile="/data2/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expGraph1Name" "$expFile" "$expGraph2Name" "$outputFile" "$relGraph1Name"
