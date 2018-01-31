# SECTION 1: recompile analysis code

make -j
if [ "$?" != 0 ]
then
    # Make failed: abort analysis
    echo "Compilation failed - correct errors in source code."
    exit $?
fi

################################################################################

# SECTION 2: define analysis workflow

while getopts "scre" opt; do
    case ${opt} in
        s)
            runSubrunPath=true
            ;;
        c)
            runChunk=true
            ;;
        r)
            runList=true
            ;;
        e)
            eachSubrun=true
            ;;
        \?)
            # Flags unrecognized - exit script and give the user a help message
            printf "\nInvalid flag given.\n\nValid flags are:\n"
            printf "    -s (analyze a single file, given as runNumber subRunNumber)\n"
            printf "    -c (analyze a chunk of subruns, given as runNumber, first subRunNumber, last subRunNumber)\n"
            printf "    -r (analyze runs listed in ../<experiment>/runsToSort.txt)\n"
            printf "    -e (analyze each subrun separately to produce its own cross section)\n"
            exit
            ;;
    esac
done


read -r experiment<experiment.txt # find out which experimental dataset to analyze

if [ "$runList" = true ]
then
    printf "\nAnalyzing all runs in runList.txt...\n"

    detectorName=$2 # which detector should be used for cross section calculation
    if [ "$experiment" = "nickel" ]
    then
        analysisDirectoryName="/data1/analysis"
    fi

    if [ "$experiment" = "tin2" ]
    then
        analysisDirectoryName="/data2/analysis"
    fi

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
    if [ "$skip" = true ]
    then
        continue
    fi

    if [ "$eachSubrun" = true ]
    then
        ./bin/eachSubrun "$analysisDirectoryName" "$experiment" "$detectorName"
    else
        ./bin/sumAll "$analysisDirectoryName" "$experiment" "$detectorName"
    fi
fi

if [ "$runChunk" = true ]
then
    runNumber=$2
    lowSubrun=$3
    highSubrun=$4
    detectorName=$5  # which detector should be used for cross section calculation

    # Analyze a chunk of subruns in a run, specified by the run, the first subrun to
    # analyze, and the last subrun to analyze
    printf "\nAnalyzing chunk of subruns in $2...\n"

    # read input filepath and output filepath
    while read l
    do
        filepaths=($l)
        if [[ ${filepaths[0]} -le $runNumber && ${filepaths[1]} -ge $runNumber ]]
        then
            analysisDirectoryName=${filepaths[3]}
            break
        fi
    done < ../"$experiment"/filepaths.txt

    ./bin/sumChunk "$analysisDirectoryName" "$experiment" "$runNumber" "$lowSubrun" "$highSubrun" "$detectorName"

fi
