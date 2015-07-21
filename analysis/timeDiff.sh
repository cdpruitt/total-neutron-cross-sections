#!/bin/bash
# calculates differences between times from a file with one decimal time per line 
while read line
do
    if [ -n "$linePrev" ]
    then
       echo $((line - linePrev)) >> timeDiff.txt
    fi
    linePrev=$line

done < timeTags.txt
