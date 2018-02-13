#!/bin/bash

expFile="/data2/analysis/total.root"
expCGraphName="CNat"
expSn112GraphName="Sn112"
expSnGraphName="SnNat"
expSn124GraphName="Sn124"
expPbGraphName="PbNat"

litFile="/data2/analysis/literatureData.root"
litCGraphName="Natural C (n,tot)"
litSn112GraphName="Sn112 (n,tot)"
litSnGraphName="Natural Sn (n,tot)"
litSn124GraphName="Sn124 (n,tot)"
litPbGraphName="Natural Pb (n,tot)"

relCGraphName="CNat, expLit"
relSn112GraphName="Sn112, expLit"
relSnGraphName="SnNat, expLit"
relSn124GraphName="Sn124, expLit"
relPbGraphName="PbNat, expLit"

outputFile="/data2/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expCGraphName" "$litFile" "$litCGraphName" "$outputFile" "$relCGraphName"
../bin/relativeDiffCS "$expFile" "$expSn112GraphName" "$litFile" "$litSn112GraphName" "$outputFile" "$relSn112GraphName"
../bin/relativeDiffCS "$expFile" "$expSnGraphName" "$litFile" "$litSnGraphName" "$outputFile" "$relSnGraphName"
../bin/relativeDiffCS "$expFile" "$expSn124GraphName" "$litFile" "$litSn124GraphName" "$outputFile" "$relSn124GraphName"
../bin/relativeDiffCS "$expFile" "$expPbGraphName" "$litFile" "$litPbGraphName" "$outputFile" "$relPbGraphName"

# scale to percentage
../bin/multiplyCS "$outputFile" "$relCGraphName" "100" "$relCGraphName, percent"
../bin/multiplyCS "$outputFile" "$relSn112GraphName" "100" "$relSn112GraphName, percent"
../bin/multiplyCS "$outputFile" "$relSnGraphName" "100" "$relSnGraphName, percent"
../bin/multiplyCS "$outputFile" "$relSn124GraphName" "100" "$relSn124GraphName, percent"
../bin/multiplyCS "$outputFile" "$relPbGraphName" "100" "$relPbGraphName, percent"
