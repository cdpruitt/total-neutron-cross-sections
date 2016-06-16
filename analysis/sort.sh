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

# SECTION 0: Create a blacklist for 'bad' runs

declare -a blacklist=("150-0011" "150-0012" "150-0013" "150-0014"
                      "150-0015" "150-0016" "150-0017" "150-0018"
                      "150-0019"
                      
                      "151-0000"
                      
                      "152-0003" "152-0004" "152-0005" "152-0006"
                      "152-0007" "152-0008" "152-0009" "152-0010"
                      "152-0011"
                      
                      "155-0008"
                      
                      "168-0008" "168-0009" "168-0010"
                      
                      "170-0003" "170-0021") 

################################################################################

# SECTION 1: recompile analysis code

# Check to see if analysis codes have been modified since last compile
if [ raw.cpp -nt raw ] || [ resort.cpp -nt resort ] || [ histos.cpp -nt histos ] ||
    [ waveform.cpp -nt waveform ] || [ sumRun.cpp -nt sumRun ] || [ sumAll.cpp -nt sumAll ]
then
    # sorting code modified
    # recompile to prepare for event sorting
    make
if [ raw.cpp -nt raw ] || [ resort.cpp -nt resort ] || [ histos.cpp -nt histos ] ||
    [ waveform.cpp -nt waveform ] || [ sumRun.cpp -nt sumRun ] || [ sumAll.cpp -nt sumAll ]
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
    # Check to make sure that the run to be sorted isn't on the blacklisted
    # runs list (SECTION 0).
    #echo "$runDir-$runNo  ${blacklist[0]}"
    for i in "${blacklist[@]}"
    do
        if [[ "$runDir-$runNo" = "$i" ]]
        then
            echo "Found sub-run \"$runDir-$runNo\" on blacklist; skipping..."
            return 1
        fi
    done

    # Convert selected .evt file to ROOT tree, stored as runX-YYYY_raw.root
    ./raw "$runDir" "$runNo" "$datapath" "$outpath"
    # If ./raw returned an exit status indicating failure, exit the script
    rc=$?; if [[ $rc != 0 ]]; then echo "./raw $runDir $runNo returned $rc, indicating error; exiting"; exit $rc; fi

    # Process raw ROOT tree into channel-specific subtrees in runX-YYYY_sorted.root
    ./resort "$runDir" "$runNo" "$outpath"
    # If ./resort returned an exit status indicating failure, exit the script
    rc=$?; if [[ $rc != 0 ]]; then echo "./resort $runDir $runNo returned $rc, indicating error; exiting"; exit $rc; fi

    # Process channel-specific subtrees into histograms in runX-YYYY_histos.root
    ./histos "$runDir" "$runNo" "$outpath"
    # If ./resort returned an exit status indicating failure, exit the script
    rc=$?; if [[ $rc != 0 ]]; then echo "./histos $runDir $runNo returned $rc, indicating error; exiting"; exit $rc; fi

    # Process waveforms from channel-specific subtrees
    ./waveform "$runDir" "$runNo" "$outpath"
    # If ./waveform returned an exit status indicating failure, exit the script
    rc=$?; if [[ $rc != 0 ]]; then echo "./waveform $runDir $runNo returned $rc, indicating error; exiting"; exit $rc; fi
}

prepDir ()
{
    # If the analysis directory in outpath doesn't exist yet, make it
    if [ ! -d $outpath/analysis/run$runDir ]
    then
        mkdir $outpath/analysis/run$runDir
    fi
}

################################################################################

# SECTION 3: process script flags to identify which runs should be analyzed

# Default state is that runs will NOT be read from a text file that lists them
runlist=false

# Parse runtime flags
while getopts "trsa" opt; do
    case ${opt} in
        t)
            # Have ./raw produce text output (currently disabled)
            printf "\nText output enabled (this will slow sorting considerably)\n"
            text=true
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

# Sort only a single sub-run, listed as arguments to this script
if [ "$single" = true ]
then
    runDir=$2
    runNo=$3

    # set file path to fine this run number's data
    if [ $2 -lt 6 ]
    then
        datapath=/media/cdpruitt/Drive3
        outpath=/data3
    elif [ $2 -gt 128 ] && [ $2 -lt 160 ]
    then
        datapath=/media/cdpruitt/Drive2
        outpath=/data2
    elif [ $2 -gt 159 ] && [ $2 -lt 178 ]
    then
        datapath=/media/cdpruitt/Drive3
        outpath=/data3
    else
        printf "\nRun directory outside bounds (runs 128-177) - check run number"
        exit
    fi

    # Prepare analysis directory
    prepDir

    runName="$datapath/output/run$runDir/data-$runNo.evt"
    printf "\nSorting single sub-run $runName\n"

    # Start sort
    runSize=$(du -k "$runName" | cut -f 1)
    if [ $runSize -ge 1000000 ]
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
    runDir=$(ls -t $datapath/output | head -1)
    runNo=$(ls -t $datapath/output/run$runDir/data-* | head -1 | egrep\
        -o '[0-9]+' | tail -1)
    runName="$datapath/output/run$runDir/data-$runNo.evt"
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
    while read runDir; do
        printf "Reading from directory run$runDir\n"

        # set path to data/analysis depending on run number
        if [ $runDir -lt 6 ]
        then
            datapath=/media/cdpruitt/Drive3
            outpath=/data3
        elif [ $runDir -gt 128 ] && [ $runDir -lt 160 ]
        then
            datapath=/media/cdpruitt/Drive2
            outpath=/data2
        elif [ $runDir -gt 159 ] && [ $runDir -lt 178 ]
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
            for f in $datapath/output/run$runDir/data-*;
            do
                runNo=$(echo $f | egrep -o '[0-9]+' | tail -1)
                runName="$datapath/output/run$runDir/data-$runNo.evt"
                printf "\n***Starting sort of sub-run $runNo***\n"

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
            ./sumRun $runDir $outpath
        else
            # Sort just the most recent sub-run in the specified run directory
            runNo=$(ls -t $datapath/output/run$runDir/data-* | head -1 | egrep\
                -o '[0-9]+' | tail -1)
            runName="$datapath/output/run$runDir/data-$runNo.evt"
            printf "\n***Starting sort of sub-run $runNo***\n"

            # Sort sub-run
            sort
        fi
    done < runsToSort.txt
fi
