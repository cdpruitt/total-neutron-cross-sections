#!/bin/bash

scaledown=$1
percent=$2

expFile="/data2/analysis/total.root"
expGraph1Name="ONat"
expGraph2Name="O18"

litFile="/data2/analysis/literatureData.root"
litGraph1Name="NatO(n,tot)"
litGraph2Name="18O(n,tot)"

relGraph1Name="Natural O, expLit"
relGraph2Name="O18, expLit"

outputFile="/data2/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expGraph1Name" "$litFile" "$litGraph1Name" "$outputFile" "$relGraph1Name"
../bin/relativeDiffCS "$expFile" "$expGraph2Name" "$litFile" "$litGraph2Name" "$outputFile" "$relGraph2Name"
