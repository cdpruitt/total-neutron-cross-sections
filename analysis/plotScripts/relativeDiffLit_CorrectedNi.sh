#!/bin/bash

expFile="/data1/analysis/corrected.root"
expCGraphName="CNat"
expNi58GraphName="Ni58"
expNiGraphName="NiNat"
expNi64GraphName="Ni64"
expPbGraphName="PbNat"

litFile="/data1/analysis/literatureData.root"
litCGraphName="Natural C (n,tot)"
litNi58GraphName="Ni58 (n,tot)"
litNiGraphName="Natural Ni (n,tot)"
litNi64GraphName="Ni64 (n,tot)"
litPbGraphName="Natural Pb (n,tot)"

relCGraphName="CNat, expLit, corrected"
relNi58GraphName="Ni58, expLit, corrected"
relNiGraphName="NiNat, expLit, corrected"
relNi64GraphName="Ni64, expLit, corrected"
relPbGraphName="PbNat, expLit, corrected"

outputFile="/data1/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expCGraphName" "$litFile" "$litCGraphName" "$outputFile" "$relCGraphName"
../bin/relativeDiffCS "$expFile" "$expNi58GraphName" "$litFile" "$litNi58GraphName" "$outputFile" "$relNi58GraphName"
../bin/relativeDiffCS "$expFile" "$expNiGraphName" "$litFile" "$litNiGraphName" "$outputFile" "$relNiGraphName"
../bin/relativeDiffCS "$expFile" "$expNi64GraphName" "$litFile" "$litNi64GraphName" "$outputFile" "$relNi64GraphName"
../bin/relativeDiffCS "$expFile" "$expPbGraphName" "$litFile" "$litPbGraphName" "$outputFile" "$relPbGraphName"

# scale to percentage
../bin/multiplyCS "$outputFile" "$relCGraphName" "100" "$relCGraphName, percent"
../bin/multiplyCS "$outputFile" "$relNi58GraphName" "100" "$relNi58GraphName, percent"
../bin/multiplyCS "$outputFile" "$relNiGraphName" "100" "$relNiGraphName, percent"
../bin/multiplyCS "$outputFile" "$relNi64GraphName" "100" "$relNi64GraphName, percent"
../bin/multiplyCS "$outputFile" "$relPbGraphName" "100" "$relPbGraphName, percent"
