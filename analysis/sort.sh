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

while getopts ":t" opt; do
    case $opt in
        t)
            printf "\nText output enabled (this will slow sorting considerably)\n"
            text=true
            ;;
        rl)
            printf "\nRunlist mode enabled. Run numbers will be read from ./runsToSort.txt"
            runlist=true
            ;;
        \?)
            printf "\nInvalid option: -$OPTARG.\n\nValid options are:"
            printf "\n    -t (text output)"
            printf "\n    -rl (use runlist)\n\n"
            exit
            ;;
    esac
done

if [ "$runlist" = false ]
then
    runNo=$(ls -t ../output | head -1 | egrep -o '[0-9]+')
    runName="../output/run$runNo.evt"
    printf "\nSorting most recent file $runName\n"
    ./sort "$runNo" "$text"
else
    if [ -a ./runsToSort.txt ]
    then
        printf "\n"
        while read runNo; do
            runName="../output/run$runNo.evt"
            printf "Sorting $runName\n"
            ./sort "$runNo" "$text"
        done < runsToSort.txt
    fi
fi
