#!/bin/bash

expFile="/data1/analysis/total.root"
expGraph1Name="Ni64PurityCorrected"
expGraph2Name="Ni58"

relGraph1Name="Ni64Ni58"

outputFile="/data1/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expGraph1Name" "$expFile" "$expGraph2Name" "$outputFile" "$relGraph1Name"

../bin/multiplyCS "$outputFile" "$relGraph1Name" "100" "$relGraph1Name, percent"
