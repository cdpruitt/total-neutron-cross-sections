#!/bin/bash

expFile="/data1/analysis/total.root"
expCGraphName="CNat"
expNi58GraphName="Ni58"
expNiGraphName="NiNat"
expNi64GraphName="Ni64PurityCorrected"
expPbGraphName="PbNat"

litFile="/data1/analysis/literatureData.root"
litCGraphName="Natural C (n,tot)"
litNi58GraphName="Ni58 (n,tot)"
litNiGraphName="Natural Ni (n,tot)"
litNi64GraphName="Ni64 (n,tot)"
litPbGraphName="Natural Pb (n,tot)"

relCGraphName="CNat, expLit, absolute"
relNi58GraphName="Ni58, expLit, absolute"
relNiGraphName="NiNat, expLit, absolute"
relNi64GraphName="Ni64, expLit, absolute"
relPbGraphName="PbNat, expLit, absolute"

relativeFile="/data1/analysis/relative.root"

../bin/relativeCS "$expFile" "$expCGraphName" "$litFile" "$litCGraphName" "$relativeFile" "$relCGraphName"
../bin/relativeCS "$expFile" "$expNi58GraphName" "$litFile" "$litNi58GraphName" "$relativeFile" "$relNi58GraphName"
../bin/relativeCS "$expFile" "$expNiGraphName" "$litFile" "$litNiGraphName" "$relativeFile" "$relNiGraphName"
../bin/relativeCS "$expFile" "$expNi64GraphName" "$litFile" "$litNi64GraphName" "$relativeFile" "$relNi64GraphName"
../bin/relativeCS "$expFile" "$expPbGraphName" "$litFile" "$litPbGraphName" "$relativeFile" "$relPbGraphName"
