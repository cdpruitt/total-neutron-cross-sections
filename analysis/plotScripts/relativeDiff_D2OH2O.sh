#!/bin/bash

expFile="/data2/analysis/total.root"
expGraph1Name="ONat_fromD2O"
expGraph2Name="ONat_fromH2O"

relGraph1Name="ONatD2O_ONatH2O"

outputFile="/data2/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expGraph1Name" "$expFile" "$expGraph2Name" "$outputFile" "$relGraph1Name"
