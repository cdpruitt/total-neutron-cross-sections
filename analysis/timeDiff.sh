#!/bin/bash
# calculates differences between times from a file with one decimal time per line 
if [ -f timeDiff.txt ]
then
    rm timeDiff.txt
fi

while read line
do
    if [ -n "$linePrev" ]
    then
        echo $((2*(line - linePrev))) >> timeDiff.txt
    fi
    linePrev=$line

done < timeTags.txt
