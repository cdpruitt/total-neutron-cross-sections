#!/bin/bash
runStem=output/run
i=0
while [[ -e $runStem$i.evt ]]
do
    let i++
done
runName=$runStem$i.evt

runMeta=$runStem$i.txt
touch $runMeta

echo "New run will be $runName."

echo "Enter run type:"
read runtype

echo "Enter run description/comments:\n"
read comments

echo "Enter channels enabled \(i.e., 0 2 4\):"
read channels

echo "Enter experimenters on shift:" 
read people

echo "RUN $i metadata\n" >> $runMeta
echo "Type: $runtype\n" >> $runMeta
echo "Description: $comments\n" >> $runMeta
echo "Channels: $channels\n" >> $runMeta
echo "On shift: $people\n" >> $runMeta

echo "Time start: $(date +%R)\n" >> $runMeta
rstart=$(date +%s)

bin/readout testConfig/config.txt $runName

echo "Time stop: $(date +%R)\n" >> $runMeta
rstop=$(date +%s)
diff=$(($rstop-$rstart))

echo "Run duration: $diff seconds" >> $runMeta
