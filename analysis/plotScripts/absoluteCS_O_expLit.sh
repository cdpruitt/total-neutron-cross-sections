#!/bin/bash

scaledown=$1

expFile="/data3/analysis/total.root"
expH2OGraphName="H2O"
expH218OGraphName="H218O"

litFile="/data3/analysis/literatureData.root"
litNatOGraphName="NatO(n,tot)"
litO18GraphName="18O(n,tot)"
litHGraphName="Hydrogen(n,tot)"

expNatOGraphName="NatO"
expO18GraphName="O18"

../subtractCS "$expFile" "$expH2OGraphName" "$litFile" "$litHGraphName" 2 1 "$expNatOGraphName"
../subtractCS "$expFile" "$expH218OGraphName" "$litFile" "$litHGraphName" 2 1 "$expO18GraphName"

if [ "$scaledown" -gt 1 ]
then
    printf "Scaling down relative CS plots by factor of $scaledown.\n"

    scaledownFile="/data3/analysis/scaledown.root"

    ../scaledownCS "$expFile" "$expH2OGraphName" "$scaledown" "$scaledownFile" "$expH2OGraphName"
    ../scaledownCS "$expFile" "$expH218OGraphName" "$scaledown" "$scaledownFile" "$expH218OGraphName"
fi
