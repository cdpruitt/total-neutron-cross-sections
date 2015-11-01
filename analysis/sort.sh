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
    runDir=$(ls -t /media/ExternalDrive1/output | head -1)
    runNo=$(ls -t /media/ExternalDrive1/output/$runDir/data-* | head -1 | egrep\
    -o '[0-9]+' | tail -1)
    runName="/media/ExternalDrive1/output/$runDir/data-$runNo.evt"
    printf "\nSorting most recent file $runName\n"

    if [ ! -a /media/ExternalDrive1/analysis/$runDir ]
    then
        mkdir /media/ExternalDrive1/analysis/$runDir
    fi
    ./sort "$runDir" "$runNo" "$text" "$cs"
else
    if [ -a ./runsToSort.txt ]
    then
        printf "\n"
        while read runDir; do
	    printf "Reading from directory run$runDir\n"

            if [ ! -a /media/ExternalDrive1/analysis/run$runDir ]
            then
                mkdir /media/ExternalDrive1/analysis/run$runDir
            fi

            if [ "$allFiles" = true ]
            then
                for f in /media/ExternalDrive1/output/run$runDir/data-*;
                do
                    runNo=$(echo $f | egrep -o '[0-9]+' | tail -1)
                    runName="/media/ExternalDrive1/output/$runDir/data-$runNo.evt"
                    printf "Sorting $runName\n"
                    ./sort "run$runDir" "$runNo" "$text" "$cs"
                done
            else
                runNo=$(ls -t /media/ExternalDrive1/output/run$runDir/data-* | head -1 | egrep\
                    -o '[0-9]+' | tail -1)
                runName="/media/ExternalDrive1/output/$runDir/data-$runNo.evt"
                printf "Sorting $runName\n"
                ./sort "run$runDir" "$runNo" "$text" "$cs"
            fi
        done < runsToSort.txt
    fi
fi
