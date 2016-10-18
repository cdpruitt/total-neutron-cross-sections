#!/bin/bash

################################################################################
#                               analyze.sh                                     #
################################################################################

# This script initiates data analysis. It decides which experimental
# runs to analyze, then calls C++ programs to process those runs.
# Several flags can be used with this script to manage this workflow:
#
#  Flag |                           Description
#-------+-----------------------------------------------------------------------
#    -f | analyze a single event file, specified by the full filepath
#       | (e.g., ./analyze -f path/to/input/file path/to/output/file)
#-------+-----------------------------------------------------------------------
#    -s | analyze a single event file, specified by run and subrun numbers
#       | (e.g., ./analyze -s <run number> <subrun number>)
#-------+-----------------------------------------------------------------------
#    -r | for each run given in ../<experiment>/runsToSort.txt, analyze the most
#       | recent subrun
#       | (e.g., ./analyze -r)
#-------+-----------------------------------------------------------------------
#    -a | if using -r, analyze ALL subruns in each run, not just the most recent
#       | (e.g., ./analyze -ra)
#-------+-----------------------------------------------------------------------
#    -t | produce a formatted text file for each data channel, instead of analysis
#       | (e.g., ./analyze -t)
#-------+-----------------------------------------------------------------------
#    -o | overwrite previous analysis histograms (use if analysis code has changed)
#       | (e.g., ./analyze -o)
#-------+-----------------------------------------------------------------------
#    -w | overwrite previous waveform fitting (use if fitting code has changed)
#       | (e.g., ./analyze -o)
#
# Flags can be combined for additional functionality (e.g., ./analyze -ro)
#
################################################################################

# SECTION 1: recompile analysis code

make
if [ "$?" != 0 ] # "$?" is return value of previous command
then
    echo "Compilation failed - correct errors in source code."
    exit $?
fi

################################################################################

# SECTION 2: process script flags

# Default analysis behavior is:
#   - to analyze only the most-recently-modified subrun
#   - not to produce text output of event data

# Parse runtime flags
while getopts "fsratow" opt; do
    case ${opt} in
        f)
            fullFilePath=true
            ;;
        s)
            runSubrunPath=true
            ;;
        r)
            runlist=true
            ;;
        a)
            allSubruns=true
            ;;
        t)
            produceText=true
            ;;
        o)  overwriteHistos=true
            ;;
        w)  overwriteWaveforms=true
            ;;
        \?)
            # Flags unrecognized - exit script and give user help text
            printf "\nInvalid option: -$OPTARG.\n\nValid options are:"
            printf "\n    -f (analyze a single file, given as filepath)\n"
            printf "\n    -s (analyze a single file, given as runNumber subRunNumber)\n"
            printf "\n    -r (analyze runs listed in ../<experiment>/runsToSort.txt)\n"
            printf "\n    -a (if using -r, read ALL subruns in each run, not just the most recent)\n"
            printf "\n    -t (produce text output of event data instead of doing full analysis)\n"
            printf "\n    -o (overwrite existing event histograms - use if analysis code has been changed)\n"
            printf "\n    -o (overwrite existing waveform fits - use if waveform fitting code has been changed)\n"

            exit
            ;;
    esac
done

################################################################################

# SECTION 3: define analysis workflow

# This function is called once for each data set we want to analyze.

analyze ()
{
    inputFileName=$1
    outputDirectoryName=$2
    runNumber=$3
    subrunNumber=$4

    if [ "$produceText" == true ]
    then
        printf "\nText output enabled... \n"
        ./text $inputFileName
        exit # we just want the text output, so stop analysis here
    fi

    # Create directory to hold the output of our analysis
    if [ ! -d $outputDirectoryName ]
    then
        mkdir $outputDirectoryName
    fi

    if [ "$overwriteHistos" == true ]
    then
        printf "\nOverwriting existing histogram file $outputDirectoryName"histos.root"...\n"
        rm $outputDirectoryName"/histos.root"
    fi

    if [ "$overwriteWaveforms" == true ]
    then
        printf "\nOverwriting existing waveform file $outputDirectoryName"waveform.root"...\n"
        rm $outputDirectoryName"/waveform.root"
    fi

    # Send error messages to a text file
    2>&1 | tee > $outputDirectoryName"/error.txt"

    ./driver $inputFileName $outputDirectoryName $runNumber

    if [ -s $outputDirectoryName"/temp.root" ]
    then
        rm $outputDirectoryName"/temp.root"
    fi
}

################################################################################

# SECTION 4: determine which runs should be analyzed

read -r experiment<experiment.txt # find out which experimental dataset to analyze

# Analyze a single event file, specified by the full filepath
if [ "$fullFilePath" = true ]
then
    printf "\nAnalyzing single event file $2...\n"
    inputFileName=$2
    outputDirectoryName=$3
    analyze $inputFileName $outputDirectoryName
    exit
fi

# Analyze a single event file, specified by its run number and subrun number
if [ "$runSubrunPath" = true ]
then
    runNumber=$2
    subrunNo=$3

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
    done < ../$experiment/filepaths.txt

    if [[ $datapath == "" ]]
    then
        echo "Failed to find filepath to input data. Exiting..."
        exit
    fi

    inputFileName="$datapath/output/run$runNumber/data-$subrunNo.evt"
    outputDirectoryName="$outpath/analysis/run$runNumber/$subrunNo/"

    printf "\nSorting single sub-run $inputFileName\n"

    # Start analysis
    analyze $inputFileName $outputDirectoryName $runNumber
    exit
fi

# Analyze runs listed in runsToSort.txt
if [[ $runlist = true && -a ../$experiment/runsToSort.txt ]]
then
    printf "\nRunlist mode enabled. Reading runs from ./runsToSort.txt..."

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
            echo $outpath
        done < ../$experiment/filepaths.txt

        if [ $datapath == 0 ]
        then
            echo "Failed to find filepath to input data. Exiting..."
            exit
        fi

        # Check to see if all subruns should be analyzed, or just one
        if [ "$allSubruns" = true ]
        then
            # Sort all sub-runs in the specified run directory
            printf "Reading all subruns in specified run\n"
            for f in $datapath/output/run$runNumber/data-*;
            do
                subrunNo=$(echo $f | egrep -o '[0-9]+' | tail -1)
                inputFileName="$datapath/output/run$runNumber/data-$subrunNo.evt"
                outputDirectoryName="$outpath/analysis/run$runNumber/$subrunNo/"
                printf "\n***Starting sort of sub-run $subrunNo***\n"

                # Skip subruns < 1GB in size
                runSize=$(du -k "$inputFileName" | cut -f 1)
                if [ $runSize -lt 1000000 ]
                then
                    printf "$inputFileName is less than 1 GB in size; ignoring...\n"
                    continue
                fi

                # Skip subruns on the blacklist
                skip=false
                while read l
                do
                    if [[ "$runNumber-$subrunNo" == "$l" ]]
                    then
                        echo "Found sub-run "$l" on blacklist; skipping..."
                        skip=true
                        break
                    fi
                done < ../$experiment/blacklist.txt
                if [ $skip == true ]
                then
                    continue
                fi

                analyze $inputFileName $outputDirectoryName $runNumber
            done

            # Last, sum together histograms from all the runs just sorted
            ./sumRun $runNumber $outpath
        else
            # Sort just the most recent sub-run in the specified run directory
            subrunNo=$(ls -t $datapath/output/run$runNumber/data-* | head -1 | egrep\
                -o '[0-9]+' | tail -1)
            inputFileName="$datapath/output/run$runNumber/data-$subrunNo.evt"
            outputDirectoryName="$outpath/analysis/run$runNumber/$subrunNo/"
            printf "\n***Starting sort of sub-run $subrunNo***\n"

            # Sort sub-run
            analyze $inputFileName $outputDirectoryName
        fi
    done < ../$experiment/runsToSort.txt
    ./sumAll $experiment
    exit
fi

# Default behavior: sort only the most recent run (as determined by files
# modified in defaultFilepath.txt in the experiment directory)

echo $experiment
read -r filepaths<../$experiment/defaultFilepath.txt
datapath=${filepaths[0]}
outpath=${filepaths[1]}

# Find highest-numbered run directory in datapath
runNumber=$(ls -t $datapath/output | head -1)
subrunNo=$(ls -t $datapath/output/run$runNumber/data-* | head -1 | egrep\
    -o '[0-9]+' | tail -1)
inputFileName="$datapath/output/run$runNumber/data-$subrunNo.evt"
outputDirectoryName="$outpath/analysis/run$runNumber/$subrunNo/"
printf "\nSorting most recent run $outputDirectoryName\n"

# Start sort
analyze $inputFileName $outputDirectoryName
