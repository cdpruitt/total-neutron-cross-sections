#!/bin/bash

expFile="/data2/analysis/total.root"
expGraph1Name="Sn124"
expGraph2Name="Sn112"

relGraph1Name="Sn"

outputFile="/data2/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expGraph1Name" "$expFile" "$expGraph2Name" "$outputFile" "$relGraph1Name"
