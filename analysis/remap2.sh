#!/bin/bash

files="[0-9]*/*.root"
regex="([0-9]*)/([a-z]).root"

for f in $files
do
    if [[ $f =~ $regex ]]
    then
        subrun="${BASH_REMATCH[1]}"
        analysisName="${BASH_REMATCH[2]}"
        echo $subrun $analysisName
        if [[ $analysisName == "s" ]]
        then
            mv $f $subrun"/sorted.root"
        fi
    fi
done
