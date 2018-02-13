#!/bin/bash

expFile="/data2/analysis/corrected.root"
expGraph1Name="Sn124"
expGraph2Name="Sn112"

relGraph1Name="Sn124Sn112Corrected"

outputFile="/data2/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expGraph1Name" "$expFile" "$expGraph2Name" "$outputFile" "$relGraph1Name"
