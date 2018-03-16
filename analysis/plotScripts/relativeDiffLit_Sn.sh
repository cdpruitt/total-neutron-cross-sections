#!/bin/bash

expFile="/data2/analysis/total.root"

expCGraphName="CNat"
expSn112GraphName="Sn112"
expSnGraphName="SnNat"
expSn124GraphName="Sn124"
expPbGraphName="PbNat"

expCGraphSEName="CNatSysErrors"
expSn112GraphSEName="Sn112SysErrors"
expSnGraphSEName="SnNatSysErrors"
expSn124GraphSEName="Sn124SysErrors"
expPbGraphSEName="PbNatSysErrors"

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

relCGraphSEName="CNat, expLitSysErrors"
relSn112GraphSEName="Sn112, expLitSysErrors"
relSnGraphSEName="SnNat, expLitSysErrors"
relSn124GraphSEName="Sn124, expLitSysErrors"
relPbGraphSEName="PbNat, expLitSysErrors"

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

../bin/relativeDiffCS "$expFile" "$expCGraphSEName" "$litFile" "$litCGraphName" "$outputFile" "$relCGraphSEName"
../bin/relativeDiffCS "$expFile" "$expSn112GraphSEName" "$litFile" "$litSn112GraphName" "$outputFile" "$relSn112GraphSEName"
../bin/relativeDiffCS "$expFile" "$expSnGraphSEName" "$litFile" "$litSnGraphName" "$outputFile" "$relSnGraphSEName"
../bin/relativeDiffCS "$expFile" "$expSn124GraphSEName" "$litFile" "$litSn124GraphName" "$outputFile" "$relSn124GraphSEName"
../bin/relativeDiffCS "$expFile" "$expPbGraphSEName" "$litFile" "$litPbGraphName" "$outputFile" "$relPbGraphSEName"

# scale to percentage
../bin/multiplyCS "$outputFile" "$relCGraphSEName" "100" "$relCGraphSEName, percent"
../bin/multiplyCS "$outputFile" "$relSn112GraphSEName" "100" "$relSn112GraphSEName, percent"
../bin/multiplyCS "$outputFile" "$relSnGraphSEName" "100" "$relSnGraphSEName, percent"
../bin/multiplyCS "$outputFile" "$relSn124GraphSEName" "100" "$relSn124GraphSEName, percent"
../bin/multiplyCS "$outputFile" "$relPbGraphSEName" "100" "$relPbGraphSEName, percent"
