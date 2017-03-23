#!/bin/bash

files="*.root"
regex="run[0-9]*-([0-9]*)_([a-z]*).root"

for f in $files
do
    if [[ $f =~ $regex ]]
    then
        echo "$f"
        subrun="${BASH_REMATCH[1]}"
        analysisName="${BASH_REMATCH[2]}"
        mkdir $subrun
        mv $f $subrun/$analysisName".root"
    fi
done

rm run[0-9]*-[0-9]*_cross-sections.root
rm run[0-9]*-[0-9]*_error.log
