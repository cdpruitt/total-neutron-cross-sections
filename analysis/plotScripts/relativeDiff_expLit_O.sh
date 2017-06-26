#!/bin/bash

scaledown=$1
percent=$2

expFile="/data3/analysis/total.root"
expNatOGraphName="NatO"
exp18OGraphName="O18"

litFile="/data3/analysis/literatureData.root"
litNatOGraphName="NatO(n,tot)"
lit18OGraphName="18O(n,tot)"

relNatOGraphName="relDiff_natO_expLit"
rel18OGraphName="relDiff_18O_expLit"

../relativeCS "$expFile" "$expNatOGraphName" "$litFile" "$litNatOGraphName" "$relNatOGraphName"
../relativeCS "$expFile" "$exp18OGraphName" "$litFile" "$lit18OGraphName" "$rel18OGraphName"

if [ "$scaledown" -gt 1 ]
then
    printf "Scaling down relative CS plots by factor of $scaledown.\n"

    scaledownFile="/data3/analysis/scaledown.root"
    ../scaledownCS "$expFile" "$relNatOGraphName" "$scaledown" "$scaledownFile" "$relNatOGraphName"
    ../scaledownCS "$expFile" "$rel18OGraphName" "$scaledown" "$scaledownFile" "$rel18OGraphName"
fi

if [ "$percent" == true ]
then
    printf "Converting relative difference fractions to percents.\n"

    percentNatOGraphName="relDiff_natO_expLit_percent"
    percent18OGraphName="relDiff_18O_expLit_percent"

    inputFile="/data3/analysis/scaledown.root"

    ../multiplyCS "$inputFile" "$relNatOGraphName" "100" "$percentNatOGraphName"
    ../multiplyCS "$inputFile" "$rel18OGraphName" "100" "$percent18OGraphName"
fi
