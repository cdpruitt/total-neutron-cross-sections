#!/bin/bash

expFile="/data2/analysis/corrected.root"
expCGraphName="CNat"
expNiGraphName="SnNat"
expPbGraphName="PbNat"

litFile="/data2/analysis/literatureData.root"
litCGraphName="Natural C (n,tot)"
litNiGraphName="Natural Sn (n,tot)"
litPbGraphName="Natural Pb (n,tot)"

relCGraphName="NatC, expLit, corrected"
relNiGraphName="NatSn, expLit, corrected"
relPbGraphName="NatPb, expLit, corrected"

outputFile="/data2/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expCGraphName" "$litFile" "$litCGraphName" "$outputFile" "$relCGraphName"
../bin/relativeDiffCS "$expFile" "$expNiGraphName" "$litFile" "$litNiGraphName" "$outputFile" "$relNiGraphName"
../bin/relativeDiffCS "$expFile" "$expPbGraphName" "$litFile" "$litPbGraphName" "$outputFile" "$relPbGraphName"

# scale to percentage
../bin/multiplyCS "$outputFile" "$relCGraphName" "100" "$relCGraphName, percent"
../bin/multiplyCS "$outputFile" "$relNiGraphName" "100" "$relNiGraphName, percent"
../bin/multiplyCS "$outputFile" "$relPbGraphName" "100" "$relPbGraphName, percent"
