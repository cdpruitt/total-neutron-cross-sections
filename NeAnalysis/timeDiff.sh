#!/bin/bash
# calculates differences between times from a file with one decimal time per line 
if [ -f sorted/timeDiff.txt ]
then
    rm sorted/timeDiff.txt
fi

while read line
do
    if [ -n "$linePrev" ]
    then
        echo $((2*(line - linePrev))) >> sorted/timeDiff.txt
    fi
    linePrev=$line

done < sorted/timeTags.txt