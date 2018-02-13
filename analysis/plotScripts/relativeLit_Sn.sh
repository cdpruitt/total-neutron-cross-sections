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

relCGraphName="CNat, expLit, absolute"
relSn112GraphName="Sn112, expLit, absolute"
relSnGraphName="SnNat, expLit, absolute"
relSn124GraphName="Sn124, expLit, absolute"
relPbGraphName="PbNat, expLit, absolute"

relativeFile="/data2/analysis/relative.root"

../bin/relativeCS "$expFile" "$expCGraphName" "$litFile" "$litCGraphName" "$relativeFile" "$relCGraphName"
../bin/relativeCS "$expFile" "$expSn112GraphName" "$litFile" "$litSn112GraphName" "$relativeFile" "$relSn112GraphName"
../bin/relativeCS "$expFile" "$expSnGraphName" "$litFile" "$litSnGraphName" "$relativeFile" "$relSnGraphName"
../bin/relativeCS "$expFile" "$expSn124GraphName" "$litFile" "$litSn124GraphName" "$relativeFile" "$relSn124GraphName"
../bin/relativeCS "$expFile" "$expPbGraphName" "$litFile" "$litPbGraphName" "$relativeFile" "$relPbGraphName"
