#!/bin/bash
runName=output/run
i=0
while [[ -e $runName$i.evt ]]
do
    let i++
done
runName=$runName$i.evt

echo "Starting "$runName"..."

bin/readout testConfig/config.txt $runName

echo "Testing after run"
