#!/bin/bash

expFile="/data2/analysis/total.root"
expGraph1Name="CNat"
expGraph2Name="shortCNat"

relGraph1Name="longCShortC"

outputFile="/data2/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expGraph1Name" "$expFile" "$expGraph2Name" "$outputFile" "$relGraph1Name"

../bin/multiplyCS "$outputFile" "$relGraph1Name" "100" "$relGraph1Name, percent"
