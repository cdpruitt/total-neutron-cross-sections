#!/bin/bash

scaledown=$1

expFile="/data2/analysis/total.root"
expSn112GraphName="Sn112"
expSnNatGraphName="NatSn"
expSn124GraphName="Sn124"

litFile="/data2/analysis/literatureData.root"
litSnGraphName="Natural Sn (n,tot)"

scaledownSn112GraphName="Sn112"
scaledownSnNatGraphName="SnNat"
scaledownSn124GraphName="Sn124"

if [ "$scaledown" -gt 1 ]
then
    printf "Scaling down relative CS plots by factor of $scaledown.\n"

    scaledownFile="/data2/analysis/scaledown.root"
    ../scaledownCS "$expFile" "$expSn112GraphName" "$scaledown" "$scaledownFile" "$scaledownSn112GraphName"
    ../scaledownCS "$expFile" "$expSnNatGraphName" "$scaledown" "$scaledownFile" "$scaledownSnNatGraphName"
    ../scaledownCS "$expFile" "$expSn124GraphName" "$scaledown" "$scaledownFile" "$scaledownSn124GraphName"
fi
