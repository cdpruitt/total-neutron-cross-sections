#!/bin/bash

if [ raw.cpp -nt raw ] || [ resort.cpp -nt resort ] || [ histos.cpp -nt histos ] # sorting code modified
then
    # recompile to prepare for event sorting
    make
    if [ raw.cpp -nt raw ] || [ resort.cpp -nt resort ] || [ histos.cpp -nt histos ]
    then
        # compilation failed - exit
        printf "\nCompilation failed - correct errors before sorting.\n"
        exit
    fi
fi

# process script options

runlist=false

while getopts "trsa" opt; do
    case ${opt} in
        t)
            printf "\nText output enabled (this will slow sorting considerably)\n"
            text=true
            ;;
        r)
            printf "\nRunlist mode enabled. Run directories will be read from ./runsToSort.txt"
            runlist=true
            ;;
        s)
            printf "\nSorting a single run: $2 $3\n"
            single=true
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
            printf "\n    -s (sort a single run using specified filename)"
            printf "\n    -a (if using runlist, read all evt files from last run)\n"
            exit
            ;;
    esac
done

if [ "$runlist" = false ]
then
    # set path to data/analysis depending on run number
    if [ $2 -lt 6 ]
    then
        datapath=/media/cdpruitt/Drive3
        outpath=/data3
    elif [ $2 -gt 128 ] && [ $2 -lt 161 ]
    then
        datapath=/media/cdpruitt/Drive2
        outpath=/data2
    elif [ $2 -gt 160 ] && [ $2 -lt 178 ]
    then
        datapath=/media/cdpruitt/Drive3
        outpath=/data3
    else
        printf "\nRun directory outside bounds (runs 128-177) - check run number"
    fi

    if [ "$single" = true ]
    then
        runDir=$2
        runNo=$3
        runName="$datapath/output/run$runDir/data-$runNo.evt"
        printf "\nSorting single file $runName\n"
    else
        runDir=$(ls -t $datapath/output | head -1)
        runNo=$(ls -t $datapath/output/run$runDir/data-* | head -1 | egrep\
            -o '[0-9]+' | tail -1)
        runName="$datapath/output/run$runDir/data-$runNo.evt"
        printf "\nSorting most recent file $runName\n"
    fi

    if [ ! -a $outpath/analysis/run$runDir ]
    then
        mkdir $outpath/analysis/run$runDir
    fi
    ./raw "$runDir" "$runNo" "$datapath" "$outpath"
    ./resort "$runDir" "$runNo" "$outpath"
    ./histos "$runDir" "$runNo" "$outpath"
else
    if [ -a ./runsToSort.txt ]
    then
        printf "\n"
        while read runDir; do
	    printf "Reading from directory run$runDir\n"

            # set path to data/analysis depending on run number
            if [ $runDir -lt 6 ]
            then
                datapath=/media/cdpruitt/Drive3
                outpath=/data3
            elif [ $runDir -gt 128 ] && [ $runDir -lt 161 ]
            then
                datapath=/media/cdpruitt/Drive2
                outpath=/data2
            elif [ $runDir -gt 160 ] && [ $runDir -lt 178 ]
            then
                datapath=/media/cdpruitt/Drive3
                outpath=/data3
            else
                printf "\nRun directory outside bounds (runs 128-177) - check run number"
            fi

            if [ ! -a $outpath/analysis/run$runDir ]
            then
                mkdir $outpath/analysis/run$runDir
            fi

            if [ "$allFiles" = true ]
            then
                for f in $datapath/output/run$runDir/data-*;
                do
                    runNo=$(echo $f | egrep -o '[0-9]+' | tail -1)
                    runName="$datapath/output/run$runDir/data-$runNo.evt"
                    printf "Sorting $runName\n"
                    ./raw "$runDir" "$runNo" "$datapath" "$outpath"
                    ./resort "$runDir" "$runNo" "$outpath"
                    ./histos "$runDir" "$runNo" "$outpath"
                done
            else
                runNo=$(ls -t $datapath/output/run$runDir/data-* | head -1 | egrep\
                    -o '[0-9]+' | tail -1)
                runName="$datapath/output/run$runDir/data-$runNo.evt"
                printf "Sorting $runName\n"
                ./raw "$runDir" "$runNo" "$datapath" "$outpath"
                ./resort "$runDir" "$runNo" "$outpath"
                ./histos "$runDir" "$runNo" "$outpath"
            fi
        done < runsToSort.txt
    fi
fi
