#!/bin/bash

litFile="/data1/analysis/literatureData.root"
litGraphName="Natural C (n,tot)"

expFile=

for i in {1..40}
do
    expFile="/data1/analysis/subRunCS/15_$i"".root"
    echo $expFile

    inputGraphName="CNat"
    outputGraphName="$inputGraphName""LitDifference"

    ../bin/subtractCS "$expFile" "$inputGraphName" "$litFile" "$litGraphName" 1 1 "$outputGraphName"
done
