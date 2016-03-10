#!/bin/bash
# converts a file with one hexidecimal time per line to a file with one decimal time per line
while read line
do
    htd $line 
done < timeTags.txt
