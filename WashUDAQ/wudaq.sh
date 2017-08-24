################################################################################
#                                  wudaq.sh                                    #
################################################################################

# This script is used to start data acquisition using the digitizer. It requires
# an installed version of the proprietary WashUDAQ software (written by Ron Fox)
# and a CAEN digitizer connected to the acquisition system.
#
# Two data acquisition modes are available:
#
#    ./wudaq.sh -p
#
#          OR
#
#    ./wudaq.sh -t
#
# ------------------------------------------------------------------------------
# | Flag  | Behavior                                                           | 
# --------+---------------------------------------------------------------------
# | -p    | Production mode. This will prompt the user for important           |
# |       | experimental details about the run (e.g., voltages of the detectors|
# |       | how many sub-runs to attempt) and save them along with the config  |
# |       | files for the run for use in analysis.                             |
# |       | Data is recorded as several 'sub-runs' (individual files) under an |
# |       | umbrella 'run' (a directory).                                      |
# --------+---------------------------------------------------------------------
# | -t    | Testing mode. This attempts to collect a single run and prompts the|
# |       | user for a location where this test run should be stored.          |
# ------------------------------------------------------------------------------
#
# For more information, visit the README in this directory.
#
################################################################################


#!/bin/bash

# Read flags to determine running mode
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

# Production mode
if [ "$production" == true ]
then
    # Make sure the staging/storage areas are accessible and have enough free
    # space
    read -r stagingLocation<"stagingLocation.txt"
    read -r storageLocation<"storageLocation.txt"
    read -r minFreeSpace<"minFreeSpace.txt"

    i=0
    while [[ -d $storageLocation$i ]]
    do
        let i++
    done

    mkdir $stagingLocation$i/

    # Ask the user for important experimental details for storing with run data
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

    # Record the config files used for this run
    cp prodConfig/config.txt $stagingLocation$i/config.txt 
    cp prodConfig/dppconfig.txt $stagingLocation$i/dppconfig.txt 
    cp prodConfig/waveformconfig.txt $stagingLocation$i/waveformconfig.txt 

    runMeta=$stagingLocation$i/meta.txt
    touch $runMeta

    echo "New run will be stored in $stagingLocation$i."

    echo "RUN $i metadata" >> $runMeta
    echo "Description: $comments" >> $runMeta
    echo "Target order: $targetOrder" >> $runMeta
    echo "Detector voltage: $detVolt" >> $runMeta
    echo "Monitor voltage: $monVolt" >> $runMeta
    echo "On shift: $people" >> $runMeta

    echo "Time start: $(date +%c)" >> $runMeta
    rstart=$(date +%s)

    runName=$stagingLocation$i/data.evt

    # run a batch-mode run (i.e., multiple subruns in a single run directory)
    # syntax is:
    # batchfile        configFile            stagingLocation         firstRunNo RunsToAttempt minFreeSpace(MB)
    bin/batchread.bash prodConfig/config.txt $stagingLocation$i/data 0          $runs         $minFreeSpace

    echo "Time stop:  $(date +%c)" >> $runMeta
    rstop=$(date +%s)

    # compute run duration
    diff=$(($rstop-$rstart))
    echo "Run duration: $diff seconds" >> $runMeta

    # list the descriptions of each run in a single file
    echo "run $i: $comments " >> "runDescriptions.log"

    # move recorded data to external disk
    mv $stagingLocation$i/ $storageLocation$i/

elif [ "$testing" == "true" ] 
    # run a single test run (not batch mode, as in production)
then

    echo "Enter location for test file:"
    read runName
    bin/readout prodConfig/config.txt $runName

else
    echo "Error: please enter a flag for the script indicating how data should 
    be collected (read documentation in wudaq.sh)"
    exit
fi
