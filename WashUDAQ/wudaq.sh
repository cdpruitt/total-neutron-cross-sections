#!/bin/bash

# describe how the run's data should be interpreted
# (e.g., 'test', 'production')
echo "Enter run type:"
read runtype

echo "Enter run description/comments:"
read comments

echo "Enter channels enabled (i.e., 0 2 4):"
read channels

echo "Enter stopping detector voltage:"
read detVolt

echo "Enter monitor detector voltage:"
read monVolt

echo "Enter experimenters on shift:" 
read people

echo "Enter number of runs to attempt:"
read runs

runStem=/media/InternalDrive2/output/run
runRepo=/media/ExternalDrive1/output/run

i=1
while [[ -d $runRepo$i ]]
do
    let i++
done

mkdir $runStem$i/

cp prodConfig/config.txt $runStem$i/config.txt 
cp prodConfig/dppconfig.txt $runStem$i/dppconfig.txt 
cp prodConfig/waveformconfig.txt $runStem$i/waveformconfig.txt 

runMeta=$runStem$i/meta.txt
touch $runMeta

echo "New run will be $runStem$i."

echo "RUN $i metadata" >> $runMeta
echo "Type: $runtype" >> $runMeta
echo "Description: $comments" >> $runMeta
echo "Channels: $channels" >> $runMeta
echo "Detector voltage: $detVolt" >> $runMeta
echo "Monitor voltage: $monVolt" >> $runMeta
echo "On shift: $people" >> $runMeta

echo "Time start: $(date +%c)" >> $runMeta
rstart=$(date +%s)

runName=$runStem$i/data.evt

# run Ron's batch run mode
# syntax is:
# batchfile configFile runStem firstRunNo RunsToAttempt minFreeSpace(MB)
bin/batchread.bash prodConfig/config.txt $runStem$i/data 0 $runs 20000 

# run a one-shot run
#bin/readout prodConfig/config.txt $runName

echo "Time stop:  $(date +%c)" >> $runMeta
rstop=$(date +%s)
diff=$(($rstop-$rstart))

echo "Run duration: $diff seconds" >> $runMeta

mv $runStem$i/ $runRepo$i/
