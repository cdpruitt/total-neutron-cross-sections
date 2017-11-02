#!/bin/bash

scaledown=$1
percent=$2

expFile="/data1/analysis/total.root"
expCGraphName="CNat"
expNiGraphName="NiNat"
expPbGraphName="PbNat"

litFile="/data1/analysis/literatureData.root"
litCGraphName="Natural carbon (n,tot)"
litNiGraphName="Natural Ni (n,tot)"
litPbGraphName="Natural Pb (n,tot)"

relCGraphName="relDiff_natC_expLit"
relNiGraphName="relDiff_natNi_expLit"
relPbGraphName="relDiff_natPb_expLit"

../bin/relativeCS "$expFile" "$expCGraphName" "$litFile" "$litCGraphName" "$relCGraphName"
../bin/relativeCS "$expFile" "$expNiGraphName" "$litFile" "$litNiGraphName" "$relNiGraphName"
../bin/relativeCS "$expFile" "$expPbGraphName" "$litFile" "$litPbGraphName" "$relPbGraphName"

if [ "$scaledown" -gt 1 ]
then
    printf "Scaling down relative CS plots by factor of $scaledown.\n"

    scaledownFile="/data1/analysis/scaledown.root"
    ../bin/scaledownCS "$expFile" "$relCGraphName" "$scaledown" "$scaledownFile" "$relCGraphName"
    ../bin/scaledownCS "$expFile" "$relNiGraphName" "$scaledown" "$scaledownFile" "$relNiGraphName"
    ../bin/scaledownCS "$expFile" "$relPbGraphName" "$scaledown" "$scaledownFile" "$relPbGraphName"
fi

if [ "$percent" == true ]
then
    printf "Converting relative difference fractions to percents.\n"

    percentCGraphName="relDiff_natC_expLit_percent"
    percentNiGraphName="relDiff_natNi_expLit_percent"
    percentPbGraphName="relDiff_natPb_expLit_percent"

    inputFile="/data1/analysis/scaledown.root"

    ../bin/multiplyCS "$inputFile" "$relCGraphName" "100" "$percentCGraphName"
    ../bin/multiplyCS "$inputFile" "$relNiGraphName" "100" "$percentNiGraphName"
    ../bin/multiplyCS "$inputFile" "$relPbGraphName" "100" "$percentPbGraphName"
fi
