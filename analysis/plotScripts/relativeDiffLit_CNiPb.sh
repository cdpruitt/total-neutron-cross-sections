#!/bin/bash

expFile="/data1/analysis/total.root"
expCGraphName="CNat"
expNiGraphName="NiNat"
expPbGraphName="PbNat"

litFile="/data1/analysis/literatureData.root"
litCGraphName="Natural carbon (n,tot)"
litNiGraphName="Natural Ni (n,tot)"
litPbGraphName="Natural Pb (n,tot)"

relCGraphName="NatC, expLit"
relNiGraphName="NatNi, expLit"
relPbGraphName="NatPb, expLit"

outputFile="/data1/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expCGraphName" "$litFile" "$litCGraphName" "$outputFile" "$relCGraphName"
../bin/relativeDiffCS "$expFile" "$expNiGraphName" "$litFile" "$litNiGraphName" "$outputFile" "$relNiGraphName"
../bin/relativeDiffCS "$expFile" "$expPbGraphName" "$litFile" "$litPbGraphName" "$outputFile" "$relPbGraphName"
