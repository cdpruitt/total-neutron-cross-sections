#!/bin/bash

expFile="/data1/analysis/corrected.root"
expGraph1Name="Ni64"
expGraph2Name="Ni58"

relGraph1Name="Ni64Ni58Corrected"

expGraphSE1Name="Ni64SysErrors"
expGraphSE2Name="Ni58SysErrors"

relGraphSE1Name="Ni64Ni58CorrectedSysErrors"

outputFile="/data1/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expGraph1Name" "$expFile" "$expGraph2Name" "$outputFile" "$relGraph1Name"

../bin/relativeDiffCS "$expFile" "$expGraphSE1Name" "$expFile" "$expGraphSE2Name" "$outputFile" "$relGraphSE1Name"
