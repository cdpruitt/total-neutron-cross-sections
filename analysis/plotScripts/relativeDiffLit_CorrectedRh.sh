#!/bin/bash

expFile="/data1/analysis/corrected.root"
expRhGraphName="RhNat"

litFile="/data1/analysis/literatureData.root"
litRhGraphName="Natural Rh (n,tot)"

relRhGraphName="RhNat, expLit, corrected"

outputFile="/data1/analysis/relative.root"

../bin/relativeDiffCS "$expFile" "$expRhGraphName" "$litFile" "$litRhGraphName" "$outputFile" "$relRhGraphName"

# scale to percentage
../bin/multiplyCS "$outputFile" "$relRhGraphName" "100" "$relRhGraphName, percent"
