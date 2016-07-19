#!/bin/bash

################################################################################
#                               sort.sh                                        #
################################################################################

# /sort.sh is used to initiate data analysis, using raw data from acquisition
# Several options can be triggered using flags when the script is run,
# modifying which runs are to be sorted. See flags section below for more info.

# There are four sections to the sort.sh script:

#   0. A list of runs to be excluded from analysis.
#   1. Check to see if the analysis code has been modified, and recompile if so
#   2. Define the analysis workflow (./raw -> ./resort -> ./histos)
#   3. Process this script's flags to determine which/how many runs to analyze
#   4. Process the desired runs/sub-runs

# Upon failure of a command, the script will be terminated and the reason for
# failure will be printed to terminal.

################################################################################

################################################################################

# SECTION 1: recompile analysis code

# Check to see if analysis codes have been modified since last compile
if [ raw.cpp -nt raw ] || [ resort.cpp -nt resort ] || [ histos.cpp -nt histos ]
    #|| [ waveform.cpp -nt waveform ] || [ sumRun.cpp -nt sumRun ]
then
    # sorting code modified
    # recompile to prepare for event sorting
    make
if [ raw.cpp -nt raw ] || [ resort.cpp -nt resort ] || [ histos.cpp -nt histos ]
    #|| [ waveform.cpp -nt waveform ] || [ sumRun.cpp -nt sumRun ]
    then
        # compilation failed - exit
        printf "\nCompilation failed - correct errors before sorting.\n"
        exit
    fi
fi

################################################################################

# SECTION 2: define analysis workflow

# This method defines the analysis workflow as raw->resort->histos.
# Once the correct runs to sort have been identified, it will be called to
# perform analysis on those runs.
sort ()
{
    # Check to make sure that the run to be sorted isn't blacklisted
    if [ -s blacklist.txt ]
    then
        while read l
        do
            if [[ "$runNumber-$subrunNo" = "$l" ]]
            then
                echo "Found sub-run "$l" on blacklist; skipping..."
                return 1
            fi
        done < blacklist.txt
    fi

    # Convert selected .evt file to ROOT tree, stored as runX-YYYY_raw.root
    ./raw "$inputFileName" "$outputFileName" 
    # If ./raw returned an exit status indicating failure, exit the script
    rc=$?; if [[ $rc != 0 ]]; then echo "./raw $runNumber $subrunNo returned $rc, indicating error; exiting"; exit $rc; fi

    # Process raw ROOT tree into channel-specific subtrees in runX-YYYY_sorted.root
    ./resort "$runNumber" "$subrunNo" "$outpath"
    # If ./resort returned an exit status indicating failure, exit the script
    rc=$?; if [[ $rc != 0 ]]; then echo "./resort $runNumber $subrunNo returned $rc, indicating error; exiting"; exit $rc; fi

    # Process channel-specific subtrees into histograms in runX-YYYY_histos.root
    ./histos "$runNumber" "$subrunNo" "$outpath"
    # If ./resort returned an exit status indicating failure, exit the script
    rc=$?; if [[ $rc != 0 ]]; then echo "./histos $runNumber $subrunNo returned $rc, indicating error; exiting"; exit $rc; fi

    # Process waveforms from channel-specific subtrees
    if [ "$waveform" = true ]
    then
        ./waveform "$runNumber" "$subrunNo" "$outpath"
        # If ./waveform returned an exit status indicating failure, exit the script
        rc=$?; if [[ $rc != 0 ]]; then echo "./waveform $runNumber $subrunNo returned $rc, indicating error; exiting"; exit $rc; fi
    fi
}

prepDir ()
{
    # If the analysis directory in outpath doesn't exist yet, make it
    if [ ! -d $outpath/analysis/run$runNumber ]
    then
        mkdir $outpath/analysis/run$runNumber
    fi
}

################################################################################

# SECTION 3: process script flags to identify which runs should be analyzed

# Default state is that runs will NOT be read from a text file that lists them
runlist=false

# Default is not to produce text output of raw run data
produceText=0

# Default is for absolute filepath NOT to be given by used
givenFile=false

# Parse runtime flags
while getopts "trsawf" opt; do
    case ${opt} in
        t)
            # Have ./raw produce text output
            printf "\nText output enabled (this will slow sorting considerably)\n"
            produceText=1
            ;;
        r)
            # Read runs to be sorted from the text file runsToSort.txt
            printf "\nRunlist mode enabled. Run directories will be read from ./runsToSort.txt"
            runlist=true
            ;;
        s)
            # Sort only one sub-run, specified after the flag as X YYYY
            single=true
            ;;
        a)
            # When the runlist is used to identify runs to sort, this will sort
            # all sub-runs in the specified run
            printf "\nReading all evt files in specified run (this will slow sorting considerably)\n"
            allFiles=true
            ;;
        w)
            # Perform a waveform fit in addition to DPP analysis
            printf "\nPerforming waveform fit in addition to DPP analysis.\n"
            waveform=true
            ;;
        f)
            # Complete path to event file given as argument
            printf "\nReading single file given by user...\n"
            givenFile=true
            ;;

        \?)
            # Flags unrecognized - exit script and give user help text
            printf "\nInvalid option: -$OPTARG.\n\nValid options are:"
            printf "\n    -t (text output)"
            printf "\n    -r (use runlist)"
            printf "\n    -s (sort a single run using specified filename)"
            printf "\n    -a (if using runlist, read all evt files from last run)\n"
            exit
            ;;
    esac
done

################################################################################

# SECTION 4: perform analysis on specified runs

# Sort a run with a specific filename
if [ "$givenFile" = true ]
then
    inputFileName=$2
    outputFileName=$3
    sort

# Sort only a single sub-run, listed as arguments to this script
elif [ "$single" = true ]
then
    runNumber=$2
    subrunNo=$3

    # set file path to fine this run number's data
    if [ "$2" -lt 6 ]
    then
        datapath=/media/cdpruitt/Drive3
        outpath=/data3
    elif [ "$2" -gt 127 ] && [ "$2" -lt 160 ]
    then
        datapath=/media/cdpruitt/Drive2
        outpath=/data2
    elif [ "$2" -gt 159 ] && [ "$2" -lt 178 ]
    then
        datapath=/media/cdpruitt/Drive3
        outpath=/data3
    else
        printf "\nRun directory outside bounds (runs 128-177) - check run number"
        exit
    fi

    # Prepare analysis directory
    prepDir

    runName="$datapath/output/run$runNumber/data-$subrunNo.evt"
    printf "\nSorting single sub-run $runName\n"

    # Start sort
    runSize=$(du -k "$runName" | cut -f 1)
    if [ "$runSize" -ge 10000 ]
    then
        sort
    else
        printf "$runName is less than 1 GB in size; ignoring...\n"
    fi

# Check to see if the runlist should be used to identify runs to sort
elif [ "$runlist" = false ]
then
    # Sort only the most recent run (as determined by the run directory number
    # in datapath)

    # We arbitrarily set datapath and outpath to Drive 3 here - edit this to use
    # another drive for sorting based on most recent run number
    datapath=/media/cdpruitt/Drive3
    outpath=/data3

    # Find highest-numbered run directory in datapath
    runNumber=$(ls -t $datapath/output | head -1)
    subrunNo=$(ls -t $datapath/output/run$runNumber/data-* | head -1 | egrep\
        -o '[0-9]+' | tail -1)
    runName="$datapath/output/run$runNumber/data-$subrunNo.evt"
    printf "\nSorting most recent run $runName\n"

    # Prepare analysis directory
    prepDir

    # Start sort
    sort

    # runlist is true, so let's check to see if a runlist actually exists
elif [ -a ./runsToSort.txt ]
then
    # The runlist exists - use it to identify runs to sort
    printf "\n"
    while read runNumber; do
        printf "Reading from directory run$runNumber\n"

        # set path to data/analysis depending on run number
        if [ $runNumber -lt 6 ]
        then
            datapath=/media/cdpruitt/Drive3
            outpath=/data3
        elif [ $runNumber -gt 127 ] && [ $runNumber -lt 160 ]
        then
            datapath=/media/cdpruitt/Drive2
            outpath=/data2
        elif [ $runNumber -gt 159 ] && [ $runNumber -lt 178 ]
        then
            datapath=/media/cdpruitt/Drive3
            outpath=/data3
        else
            printf "\nRun directory outside bounds (runs 128-177) - check run number"
        fi

        # Prepare analysis directory
        prepDir

        # Check to see if all sub-runs in this run directory should be sorted,
        # or just one
        if [ "$allFiles" = true ]
        then
            # Sort all sub-runs in the specified run directory
            for f in $datapath/output/run$runNumber/data-*;
            do
                subrunNo=$(echo $f | egrep -o '[0-9]+' | tail -1)
                runName="$datapath/output/run$runNumber/data-$subrunNo.evt"
                printf "\n***Starting sort of sub-run $subrunNo***\n"

                # Sort sub-run
                runSize=$(du -k "$runName" | cut -f 1)
                if [ $runSize -ge 1000000 ]
                then
                    sort
                else
                    printf "$runName is less than 1 GB in size; ignoring...\n"
                fi
            done

            # Last, sum together histograms from all the runs just sorted
            ./sumRun $runNumber $outpath
        else
            # Sort just the most recent sub-run in the specified run directory
            subrunNo=$(ls -t $datapath/output/run$runNumber/data-* | head -1 | egrep\
                -o '[0-9]+' | tail -1)
            runName="$datapath/output/run$runNumber/data-$subrunNo.evt"
            printf "\n***Starting sort of sub-run $subrunNo***\n"

            # Sort sub-run
            sort
        fi
    done < runsToSort.txt
fi
