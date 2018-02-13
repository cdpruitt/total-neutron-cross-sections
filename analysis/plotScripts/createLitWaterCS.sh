#!/bin/bash

litFile="/data2/analysis/literatureData.root"

litHGraphName="Hydrogen (n,tot)"
litDGraphName="Deuterium (n,tot)"

litOGraphName="NatO(n,tot)"
litO18GraphName="18O(n,tot)"

litH2OGraphName="H2O (n,tot)"
litD2OGraphName="D2O (n,tot)"
litH2O18GraphName="H2O18 (n,tot)"

../bin/subtractCS "$litFile" "$litOGraphName" "$litFile" "$litHGraphName"\
    -2 1 "$litH2OGraphName"

../bin/subtractCS "$litFile" "$litOGraphName" "$litFile" "$litDGraphName"\
    -2 1 "$litD2OGraphName"

../bin/subtractCS "$litFile" "$litO18GraphName" "$litFile" "$litHGraphName"\
    -2 1 "$litH2O18GraphName"
