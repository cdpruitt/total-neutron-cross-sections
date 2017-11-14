#!/bin/bash

scaledown=$1

expFile="/data1/analysis/relative.root"
expGraphName="relativeNi64Ni58"

if [ "$scaledown" -gt 1 ]
then
    printf "Scaling down relative CS plots by factor of $scaledown.\n"

    scaledownFile="/data1/analysis/scaledown.root"
    scaledownGraphName="relativeNi64Ni58_scaledown_$scaledown"
    ../bin/scaledownCS "$expFile" "$expGraphName" "$scaledown" "$scaledownFile" "$scaledownGraphName"
fi

printf "Converting relative difference fractions to percents.\n"

percentGraphName="relativeNi64Ni58_percent"

../bin/multiplyCS "$scaledownFile" "$scaledownGraphName" "100" "$perfectGraphName"
