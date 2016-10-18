#!/bin/bash

while getopts "pt" opt; do
    case ${opt} in
        p)
            production=true
            ;;
        t)
            testing=true
            ;;
        \?)
            # Flag unrecognized - exit script and give user help text
            printf "\nInvalid option: -$OPTARG.\n\nValid options are:"
            printf "\n    -t (take a single test run)\n"
            printf "\n    -p (take production mode data)\n"
            exit
            ;;
    esac
done

if [ "$production" == true ]
then
    # Production mode: get important experimental details

    read -r stagingLocation<"stagingLocation.txt"
    read -r storageLocation<"storageLocation.txt"

    echo "Enter run description/comments:"
    read comments

    echo "Enter target order:"
    read targetOrder

    echo "Enter stopping detector voltage:"
    read detVolt

    echo "Enter monitor detector voltage:"
    read monVolt

    echo "Enter experimenters on shift:" 
    read people

    echo "Enter number of runs to attempt:"
    read runs

    # Make staging directory on internal disk for this run
    i=0
    while [[ -d $storageLocation$i ]]
    do
        let i++
    done

    mkdir $stagingLocation$i/

    # Make copies of the config files used for this run
    cp prodConfig/config.txt $stagingLocation$i/config.txt 
    cp prodConfig/dppconfig.txt $stagingLocation$i/dppconfig.txt 
    cp prodConfig/waveformconfig.txt $stagingLocation$i/waveformconfig.txt 

    runMeta=$stagingLocation$i/meta.txt
    touch $runMeta

    echo "New run will be $stagingLocation$i."

    echo "RUN $i metadata" >> $runMeta
    echo "Description: $comments" >> $runMeta
    echo "Target order: $targetOrder" >> $runMeta
    echo "Detector voltage: $detVolt" >> $runMeta
    echo "Monitor voltage: $monVolt" >> $runMeta
    echo "On shift: $people" >> $runMeta

    echo "Time start: $(date +%c)" >> $runMeta
    rstart=$(date +%s)

    runName=$stagingLocation$i/data.evt

  # run Ron's batch run mode
  # syntax is:
  # batchfile          configFile            stagingLocation         firstRunNo RunsToAttempt minFreeSpace(MB)
    bin/batchread.bash prodConfig/config.txt $stagingLocation$i/data 0          $runs         20000 
    echo "Time stop:  $(date +%c)" >> $runMeta
    rstop=$(date +%s)
    diff=$(($rstop-$rstart))

    echo "Run duration: $diff seconds" >> $runMeta

    # list the descriptions of each run in a single file
    echo "run $i: $comments " >> "runDescriptions.log"

    # move recorded data to external disk
    mv $stagingLocation$i/ $storageLocation$i/

elif [ "$testing" == "true" ] 
# run a one-shot run
then

    echo "Enter location for test file:"
    read runName
    bin/readout prodConfig/config.txt $runName

else
    echo "Error: please enter a flag for the script indicating how data should 
be collected (read documentation in wudaq.sh)"
    exit
fi
