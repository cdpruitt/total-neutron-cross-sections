#!/bin/bash

expFile="/data2/analysis/total.root"
expCGraphName="CNat"
expSnGraphName="SnNat"
expPbGraphName="PbNat"

litFile="/data2/analysis/literatureData.root"
litCGraphName="Natural C (n,tot)"
litSnGraphName="Natural Sn (n,tot)"
litPbGraphName="Natural Pb (n,tot)"

relCGraphName="NatC, expLit"
relSnGraphName="NatSn, expLit"
relPbGraphName="NatPb, expLit"

outputFile="/data2/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expCGraphName" "$litFile" "$litCGraphName" "$outputFile" "$relCGraphName"
../bin/relativeDiffCS "$expFile" "$expSnGraphName" "$litFile" "$litSnGraphName" "$outputFile" "$relSnGraphName"
../bin/relativeDiffCS "$expFile" "$expPbGraphName" "$litFile" "$litPbGraphName" "$outputFile" "$relPbGraphName"

# scale to percentage
../bin/multiplyCS "$outputFile" "$relCGraphName" "100" "$relCGraphName, percent"
../bin/multiplyCS "$outputFile" "$relSnGraphName" "100" "$relSnGraphName, percent"
../bin/multiplyCS "$outputFile" "$relPbGraphName" "100" "$relPbGraphName, percent"
