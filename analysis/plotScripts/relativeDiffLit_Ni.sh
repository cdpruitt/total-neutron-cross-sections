#!/bin/bash

expFile="/data1/analysis/total.root"
expCGraphName="CNat"
expNi58GraphName="Ni58"
expNiGraphName="NiNat"
expNi64GraphName="Ni64PurityCorrected"
expPbGraphName="PbNat"

expCGraphSEName="CNatSysErrors"
expNi58GraphSEName="Ni58SysErrors"
expNiGraphSEName="NiNatSysErrors"
expNi64GraphSEName="Ni64PurityCorrectedSysErrors"
expPbGraphSEName="PbNatSysErrors"

litFile="/data1/analysis/literatureData.root"
litCGraphName="Natural C (n,tot)"
litNi58GraphName="Ni58 (n,tot)"
litNiGraphName="Natural Ni (n,tot)"
litNi64GraphName="Ni64 (n,tot)"
litPbGraphName="Natural Pb (n,tot)"

relCGraphName="CNat, expLit"
relNi58GraphName="Ni58, expLit"
relNiGraphName="NiNat, expLit"
relNi64GraphName="Ni64, expLit"
relPbGraphName="PbNat, expLit"

relCGraphSEName="CNat, expLit, SysErrors"
relNi58GraphSEName="Ni58, expLit, SysErrors"
relNiGraphSEName="NiNat, expLit, SysErrors"
relNi64GraphSEName="Ni64, expLit, SysErrors"
relPbGraphSEName="PbNat, expLit, SysErrors"

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

# calculate for systematic errors graphs
../bin/relativeDiffCS "$expFile" "$expCGraphSEName" "$litFile" "$litCGraphName" "$outputFile" "$relCGraphSEName"
../bin/relativeDiffCS "$expFile" "$expNi58GraphSEName" "$litFile" "$litNi58GraphName" "$outputFile" "$relNi58GraphSEName"
../bin/relativeDiffCS "$expFile" "$expNiGraphSEName" "$litFile" "$litNiGraphName" "$outputFile" "$relNiGraphSEName"
../bin/relativeDiffCS "$expFile" "$expNi64GraphSEName" "$litFile" "$litNi64GraphName" "$outputFile" "$relNi64GraphSEName"
../bin/relativeDiffCS "$expFile" "$expPbGraphSEName" "$litFile" "$litPbGraphName" "$outputFile" "$relPbGraphSEName"

# scale to percentage
../bin/multiplyCS "$outputFile" "$relCGraphSEName" "100" "$relCGraphSEName, percent"
../bin/multiplyCS "$outputFile" "$relNi58GraphSEName" "100" "$relNi58GraphSEName, percent"
../bin/multiplyCS "$outputFile" "$relNiGraphSEName" "100" "$relNiGraphSEName, percent"
../bin/multiplyCS "$outputFile" "$relNi64GraphSEName" "100" "$relNi64GraphSEName, percent"
../bin/multiplyCS "$outputFile" "$relPbGraphSEName" "100" "$relPbGraphSEName, percent"
