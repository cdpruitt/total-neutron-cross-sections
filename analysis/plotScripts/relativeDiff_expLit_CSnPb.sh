#!/bin/bash

scaledown=$1
percent=$2

expFile="/data2/analysis/total.root"
expCGraphName="longCarbon"
expSnGraphName="NatSn"
expPbGraphName="NatPb"

litFile="/data2/analysis/literatureData.root"
litCGraphName="Natural carbon (n,tot)"
litSnGraphName="Natural Sn (n,tot)"
litPbGraphName="Natural Pb (n,tot)"

relCGraphName="relDiff_natC_expLit"
relSnGraphName="relDiff_natSn_expLit"
relPbGraphName="relDiff_natPb_expLit"

../relativeCS "$expFile" "$expCGraphName" "$litFile" "$litCGraphName" "$relCGraphName"
../relativeCS "$expFile" "$expSnGraphName" "$litFile" "$litSnGraphName" "$relSnGraphName"
../relativeCS "$expFile" "$expPbGraphName" "$litFile" "$litPbGraphName" "$relPbGraphName"

if [ "$scaledown" -gt 1 ]
then
    printf "Scaling down relative CS plots by factor of $scaledown.\n"

    scaledownFile="/data2/analysis/scaledown.root"
    ../scaledownCS "$expFile" "$relCGraphName" "$scaledown" "$scaledownFile" "$relCGraphName"
    ../scaledownCS "$expFile" "$relSnGraphName" "$scaledown" "$scaledownFile" "$relSnGraphName"
    ../scaledownCS "$expFile" "$relPbGraphName" "$scaledown" "$scaledownFile" "$relPbGraphName"
fi

if [ "$percent" == true ]
then
    printf "Converting relative difference fractions to percents.\n"

    percentCGraphName="relDiff_natC_expLit_percent"
    percentSnGraphName="relDiff_natSn_expLit_percent"
    percentPbGraphName="relDiff_natPb_expLit_percent"

    inputFile="/data2/analysis/scaledown.root"

    ../multiplyCS "$inputFile" "$relCGraphName" "100" "$percentCGraphName"
    ../multiplyCS "$inputFile" "$relSnGraphName" "100" "$percentSnGraphName"
    ../multiplyCS "$inputFile" "$relPbGraphName" "100" "$percentPbGraphName"
fi
