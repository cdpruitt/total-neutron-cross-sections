#!/bin/bash
for x in $(find . -name \*sort.root); do
    mv $x $(echo "$x" | sed 's/sort/sorted/')
done
