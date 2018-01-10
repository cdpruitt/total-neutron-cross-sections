#!/bin/bash

litFile="/data2/analysis/literatureData.root"

lit18OGraphName="18O(n,tot)"

# shift O18 cross section up by 1 barn, for readability on plots

shiftFileName="/data2/analysis/literatureData.root"
shiftGraphName="$lit18OGraphName+1Barn"

../bin/shiftCS "$litFile" "$lit18OGraphName" "1" "$shiftFileName" "$shiftGraphName"
