#!/bin/bash

################################################################################
#                               analyze.sh                                     #
################################################################################
#
# Running this script initiates data analysis. It decides which experimental
# runs to analyze, then calls executables in bin/ to process those runs.
# The flags listed below specify which runs should be analyzed and the type of
# analysis to be performed on these runs.
#
# To specify which runs should be analyzed, use the following flags:
#
#  Flag |                           Description
#-------+-----------------------------------------------------------------------
#    -f | analyze a single event file, specified by the full filepath
#       | ( ./analyze -f path/to/input/file path/to/output/file)
#-------+-----------------------------------------------------------------------
#    -s | analyze a single event file, specified by run and subrun numbers
#       | ( ./analyze -s <run number> <subrun number>)
#-------+-----------------------------------------------------------------------
#    -c | analyze a chunk of subruns, specified by run, first, and last subruns
#       | ( ./analyze -c <run number> <first subrun number> <last subrun number>)
#-------+-----------------------------------------------------------------------
#    -r | analyze each run given in ../<experiment>/runsToSort.txt
#       | ( ./analyze -r)
#
# To specify how analysis should proceed, use the following flags:
#
#  Flag |                           Description
#-------+-----------------------------------------------------------------------
#    -n | normal analysis mode: start with the raw digitizer run data (.evt
#       | file) and generate cross section graphs from DPP mode event data.
#-------+-----------------------------------------------------------------------
#    -t | diagnostic mode: produce a human-readable text dump for each digitizer
#       | channel; do not produce cross sections or ROOT trees
#-------+-----------------------------------------------------------------------
#    -i | timing testing mode: produce plots showing timing coincidence between
#       | the left and right detectors for every event
#-------+-----------------------------------------------------------------------
#    -w | waveform analysis mode: 
#       | 
#-------+-----------------------------------------------------------------------
#    -d | DPP waveform analysis mode:
#       | 
#
# Examples:
#
#   ./analyze.sh -rn
#   
#   Using data from all runs listed in runsToSort.txt, produce cross sections
#   from DPP event data.
#
#
#   ./analyze.sh -ft myExperiment/1/data-0000.evt
#
#   Using the single data file myExperiment/1/data-0000.evt, produce human-
#   readable text files for each digitizer channel for inspection. No cross
#   sections, ROOT trees, or histograms will be created.
#
#   ./analyze.sh -fi myExperiment/1/data-0000.evt
#
#   Using the single data file myExperiment/1/data-0000.evt, produce plots 
#   showing timing coincidence between the two main detector arms
#
#
################################################################################

# SECTION 1: recompile analysis code


printf "*** recompile analysis code ***\n"
make -j
if [ "$?" != 0 ]
then
    # Make failed: abort analysis
    printf "Compilation failed - correct errors in source code.\n"
    exit $?
fi
printf "*** re-compilation complete ***\n"

################################################################################

# SECTION 2: process flags

fullFilePath=false;
runSubrunPath=false;
runChunk=false;
runList=false;
timeCheckMode=false;
produceText=false;
overwriteHistos=false;
useVetoPaddle=false;
overwriteWaveform=false;

while getopts "fscrithovw" opt; do
    case ${opt} in
        f)
            fullFilePath=true
            ;;
        s)
            runSubrunPath=true
            ;;
        c)
            runChunk=true
            ;;
        r)
            runList=true
            ;;
        i)
            timeCheckMode=true
            ;;
        t)
            produceText=true
            ;;
        h)  
            overwriteHistos=true
            ;;
        o)
            oneEach=true
            ;;
        v)  
            useVetoPaddle=true
            ;;
        w)  
            overwriteWaveforms=true
            ;;
        \?)
            # Flags unrecognized - exit script and give the user a help message
            printf "\nInvalid flag given.\n\nValid flags are:\n"
            printf "    -f (analyze a single file, given as filepath)\n"
            printf "    -s (analyze a single file, given as runNumber subRunNumber)\n"
            printf "    -c (analyze a chunk of subruns, given as runNumber, first subRunNumber, last subRunNumber)\n"
            printf "    -r (analyze runs listed in ../<experiment>/runsToSort.txt)\n"
            printf "    -n (normal analysis mode: produce cross sections)\n"
            printf "    -t (produce text output of event data instead of doing full analysis)\n"
            printf "    -h (normal analysis mode, but overwrite existing histograms)\n"
            printf "    -w (normal analysis mode, but overwrite existing waveform analysis)\n"

            exit
            ;;
    esac
done

################################################################################

# SECTION 3: define analysis workflow

# This function is called once for each data set we want to analyze.

analyze ()
{
    if [ "$produceText" == true ]
    then
        printf "\nText output enabled... \n"

        inputFileName=$1

        ./bin/text "$inputFileName"
        exit # we just want the text output, so stop analysis here
    fi

    if [ "$timeCheckMode" == true ]
    then
        printf "\nFine time check enabled... \n"

        inputFileName=$1
        outputDirectoryName=$2
        experiment=$3
        runNumber=$4

        ./bin/detTimeCheck "$inputFileName" "$outputDirectoryName" "$experiment" "$runNumber"
        exit # we just want the time check output, so stop analysis here
    fi

    outputDirectoryName=$2
    experiment=$3
    runNumber=$4

    # Create directory to hold the output of our analysis
    if [ ! -d "$outputDirectoryName" ]
    then
        mkdir "$outputDirectoryName"
    fi

    if [ "$overwriteHistos" == true ]
    then
        printf "\nOverwriting existing histogram file $outputDirectoryName"histos.root"...\n"
        rm "$outputDirectoryName/histos.root"
    fi

    if [ "$overwriteWaveforms" == true ]
    then
        printf "\nOverwriting existing waveform file $outputDirectoryName"waveform.root"...\n"
        rm "$outputDirectoryName/waveform.root"
    fi

    # Send error messages to a text file

    useVetoPaddle=$5
    processWaveformEvents=$6
    processDPPWaveformEvents=$7

    ./bin/driver "$inputFileName" "$outputDirectoryName" "$experiment"\
        "$runNumber" "$useVetoPaddle" "$processWaveformEvents"\
        "$processDPPWaveformEvents" 3>&1 1>&2 2>&3 \
        | tee "$outputDirectoryName"error.txt
}

################################################################################

# SECTION 4: determine which runs should be analyzed

printf "\n***** Start analysis *****\n"

experimentFileName="experiment.txt"
if [ ! -f $experimentFileName ]
then
    printf "Error: failed to find $experimentFileName. Aborting analysis.\n"
    exit
fi
read -r experiment<experiment.txt # find out which experimental dataset to analyze

if [ ! $experiment ]
then
    printf "Error: failed to read experimental directory name from\
    $experimentFileName. Is the file empty? Aborting analysis.\n"
    exit
fi

# Analyze a single event file, specified by the full filepath
if [ "$fullFilePath" = true ]
then
    printf "Full filepath mode (-f): analyzing single event file $2...\n"
    inputFileName=$2
    outputDirectoryName=$3
    runNumber=$5
    analyze "$inputFileName" "$outputDirectoryName" "$runNumber" "$useVetoPaddle" "false" "false"
    printf "Finished analysis.\n"
    exit
fi

# Analyze a single event file, specified by its run number and subrun number
if [ "$runSubrunPath" = true ]
then
    runNumber=$2
    subrunNo=$3

    printf "Single subrun mode (-s): analyzing single subrun $2...\n"

    # read input filepath and output filepath
    while read l
    do
        filepaths=($l)
        if [[ ${filepaths[0]} -le $runNumber && ${filepaths[1]} -ge $runNumber ]]
        then
            datapath=${filepaths[2]}
            outpath=${filepaths[3]}
            break
        fi
    done < ../"$experiment"/filepaths.txt

    if [[ $datapath == "" ]]
    then
        printf "Failed to find filepath to input data. Aborting analysis.\n"
        exit
    fi

    # Create directory to hold the output of our analysis
    if [ ! -d "$outpath/$runNumber" ]
    then
        mkdir "$outpath/$runNumber"
    fi

    inputFileName="$datapath/$runNumber/data-$subrunNo.evt"
    outputDirectoryName="$outpath/$runNumber/$subrunNo/"

    # Start analysis
    analyze "$inputFileName" "$outputDirectoryName" "$experiment" "$runNumber" "$useVetoPaddle" "false" "false"

    printf "Finished analysis.\n"
    exit
fi

# Analyze a chunk of subruns in a run, specified by the run, the first subrun to
# analyze, and the last subrun to analyze
if [ "$runChunk" = true ]
then
    printf "Analyzing chunk of subruns in $2...\n"
    runNumber=$2
    lowSubrun=$3
    highSubrun=$4
    subrunNo="$lowSubrun"

    while [ $subrunNo -le $highSubrun ]
    do
        # read input filepath and output filepath
        while read l
        do
            filepaths=($l)
            if [[ ${filepaths[0]} -le $runNumber && ${filepaths[1]} -ge $runNumber ]]
            then
                datapath=${filepaths[2]}
                outpath=${filepaths[3]}
                break
            fi
        done < ../"$experiment"/filepaths.txt

        if [[ $datapath == "" ]]
        then
            printf "Failed to find filepath to input data. Exiting...\n"
            exit
        fi

        # Create directory to hold the output of our analysis
        if [ ! -d "$outpath/$runNumber" ]
        then
            mkdir "$outpath/$runNumber"
        fi

        inputFileName="$datapath/$runNumber/data-$subrunNo.evt"
        outputDirectoryName="$outpath/$runNumber/$subrunNo/"

        printf "Sorting sub-run $inputFileName\n"

        # Start analysis
        analyze "$inputFileName" "$outputDirectoryName" "$experiment" "$runNumber" "$useVetoPaddle" "false" "false"
        subrunNo=$(printf "%04d" $((10#$subrunNo+1)))
    done
    exit
fi

# Analyze runs listed in runsToSort.txt
if [[ $runList = true && -a ../$experiment/runsToSort.txt ]]
then
    printf "RunList mode enabled. Reading runs from ./runsToSort.txt...\n"

    # loop through all runs listed in runsToSort.txt
    while read runNumber; do

        # read input filepath and output filepath
        while read l
        do
            IFS=' ' read -r -a tokens <<< "$l"

            if [[ "${tokens[0]}" -le "$runNumber" && "${tokens[1]}" -ge "$runNumber" ]]
            then
                datapath=${tokens[2]}
                outpath=${tokens[3]}
                break
            fi
        done < ../"$experiment"/filepaths.txt

        if [ "$datapath" == 0 ]
        then
            printf "Failed to find filepath to input data. Exiting...\n"
            exit
        fi

        # Create directory to hold the output of our analysis
        if [ ! -d "$outpath/$runNumber" ]
        then
            mkdir "$outpath/$runNumber"
        fi

        if [[ $oneEach = true ]]
        then
            inputFileName="$datapath/$runNumber/data-0000.evt"
            outputDirectoryName="$outpath/$runNumber/0000/"
            printf "\n***Starting analysis of sub-run 0000***\n"

            analyze "$inputFileName" "$outputDirectoryName" "$experiment" "$runNumber" "$useVetoPaddle" "false" "false"
            printf "\n***Finished analysis of sub-run $subrunNo***\n"

            continue
        fi

        printf "\n************************************\n"
        printf "     Reading subruns from $runNumber\n"
        printf "************************************\n"

        # Sort all sub-runs in the specified run directory
        for f in "$datapath"/"$runNumber"/data-*;
        do
            subrunNo=$(echo $f | egrep -o '[0-9]+' | tail -1)
            inputFileName="$datapath/$runNumber/data-$subrunNo.evt"

            # Skip subruns on the blacklist
            skip=false
            while read l
            do
                if [[ "$runNumber-$subrunNo" == "$l" ]]
                then
                    printf "Found sub-run "$l" on blacklist; skipping...\n"
                    skip=true
                    break
                fi
            done < ../$experiment/blacklist.txt
            if [ $skip == true ]
            then
                continue
            fi

            outputDirectoryName="$outpath/$runNumber/$subrunNo/"
            printf "\n***Starting analysis of sub-run $subrunNo***\n"

            analyze "$inputFileName" "$outputDirectoryName" "$experiment" "$runNumber" "$useVetoPaddle" "false" "false"
            printf "\n***Finished analysis of sub-run $subrunNo***\n"

        done

    done < ../"$experiment"/runsToSort.txt

    exit
fi
