#!/bin/bash

expFile="/data1/analysis/corrected.root"
expGraph1Name="Ni64"
expGraph2Name="Ni58"

relGraph1Name="Ni64Ni58Corrected"

outputFile="/data1/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expGraph1Name" "$expFile" "$expGraph2Name" "$outputFile" "$relGraph1Name"
