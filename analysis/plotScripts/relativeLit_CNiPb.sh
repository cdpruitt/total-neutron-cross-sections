#!/bin/bash

expFile="/data1/analysis/total.root"
expCGraphName="CNat"
expNiGraphName="NiNat"
expPbGraphName="PbNat"

litFile="/data1/analysis/literatureData.root"
litCGraphName="Natural carbon (n,tot)"
litNiGraphName="Natural Ni (n,tot)"
litPbGraphName="Natural Pb (n,tot)"

relCGraphName="NatC, expLit absolute"
relNiGraphName="NatNi, expLit absolute"
relPbGraphName="NatPb, expLit absolute"

relativeFile="/data1/analysis/relative.root"

../bin/relativeCS "$expFile" "$expCGraphName" "$litFile" "$litCGraphName" "$relativeFile" "$relCGraphName"
../bin/relativeCS "$expFile" "$expNiGraphName" "$litFile" "$litNiGraphName" "$relativeFile" "$relNiGraphName"
../bin/relativeCS "$expFile" "$expPbGraphName" "$litFile" "$litPbGraphName" "$relativeFile" "$relPbGraphName"
