#!/bin/bash

expFile="/data2/analysis/total.root"
expCGraphName="CNat"
expSnGraphName="SnNat"
expPbGraphName="PbNat"

litFile="/data2/analysis/literatureData.root"
litCGraphName="Natural C (n,tot)"
litSnGraphName="Natural Sn (n,tot)"
litPbGraphName="Natural Pb (n,tot)"

relCGraphName="NatC, expLit absolute"
relSnGraphName="NatSn, expLit absolute"
relPbGraphName="NatPb, expLit absolute"

relativeFile="/data2/analysis/relative.root"

../bin/relativeCS "$expFile" "$expCGraphName" "$litFile" "$litCGraphName" "$relativeFile" "$relCGraphName"
../bin/relativeCS "$expFile" "$expSnGraphName" "$litFile" "$litSnGraphName" "$relativeFile" "$relSnGraphName"
../bin/relativeCS "$expFile" "$expPbGraphName" "$litFile" "$litPbGraphName" "$relativeFile" "$relPbGraphName"
