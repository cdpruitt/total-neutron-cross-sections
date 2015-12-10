#!/bin/bash

if [ sort.cpp -nt sort ] # sort.cpp has been edited
then
    # recompile sort.cpp to prepare for event sorting
    make
    if [ sort.cpp -nt sort ]
    then
        # compilation failed - exit
        printf "\nCompilation failed - correct errors before sorting.\n"
        exit
    fi
fi

# process script options

runlist=false
runpath=/media/ExternalDrive1/output
analysispath=/media/ExternalDrive1/analysis

while getopts "trca" opt; do
    case ${opt} in
        t)
            printf "\nText output enabled (this will slow sorting considerably)\n"
            text=true
            ;;
        r)
            printf "\nRunlist mode enabled. Run directories will be read from ./runsToSort.txt"
            runlist=true
            ;;
        c)
            printf "\nCross-section histograms enabled\n"
            cs=true
            ;;
        a)
            printf "\nReading all evt files in last run (this will slow sorting\
            considerably)\n"
            allFiles=true
            ;;
        \?)
            printf "\nInvalid option: -$OPTARG.\n\nValid options are:"
            printf "\n    -t (text output)"
            printf "\n    -r (use runlist)"
            printf "\n    -c (produce cross-section plots)"
            printf "\n    -a (if using runlist, read all evt files from last run)\n"
            exit
            ;;
    esac
done

if [ "$runlist" = false ]
then
    runDir=$(ls -t $runpath | head -1)
    runNo=$(ls -t $runpath/$runDir/data-* | head -1 | egrep\
    -o '[0-9]+' | tail -1)
    runName="$runpath/$runDir/data-$runNo.evt"
    printf "\nSorting most recent file $runName\n"

    if [ ! -a $analysispath/$runDir ]
    then
        mkdir $analysispath/$runDir
    fi
    ./sort "$runDir" "$runNo" "$text" "$cs"
else
    if [ -a ./runsToSort.txt ]
    then
        printf "\n"
        while read runDir; do
	    printf "Reading from directory run$runDir\n"

            if [ ! -a $analysispath/run$runDir ]
            then
                mkdir $analysispath/run$runDir
            fi

            if [ "$allFiles" = true ]
            then
                for f in $runpath/run$runDir/data-*;
                do
                    runNo=$(echo $f | egrep -o '[0-9]+' | tail -1)
                    runName="$runpath/$runDir/data-$runNo.evt"
                    printf "Sorting $runName\n"
                    ./sort "run$runDir" "$runNo" "$text" "$cs"
                done
            else
                runNo=$(ls -t $runpath/run$runDir/data-* | head -1 | egrep\
                    -o '[0-9]+' | tail -1)
                runName="$runpath/$runDir/data-$runNo.evt"
                printf "Sorting $runName\n"
                ./sort "run$runDir" "$runNo" "$text" "$cs"
            fi
        done < runsToSort.txt
    fi
fi
