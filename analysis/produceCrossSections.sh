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

while getopts "scr" opt; do
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
        \?)
            # Flags unrecognized - exit script and give the user a help message
            printf "\nInvalid flag given.\n\nValid flags are:\n"
            printf "    -s (analyze a single file, given as runNumber subRunNumber)\n"
            printf "    -c (analyze a chunk of subruns, given as runNumber, first subRunNumber, last subRunNumber)\n"
            printf "    -r (analyze runs listed in ../<experiment>/runsToSort.txt)\n"
            exit
            ;;
    esac
done

read -r experiment<experiment.txt # find out which experimental dataset to analyze

if [ "$runList" = true ]
then
    printf "\nAnalyzing all runs in runList.txt...\n"

    analysisDirectoryName="/data1/analysis"
    ./bin/sumAll "$analysisDirectoryName" "$experiment"
fi

if [ "$runChunk" = true ]
then
    # Analyze a chunk of subruns in a run, specified by the run, the first subrun to
    # analyze, and the last subrun to analyze
    printf "\nAnalyzing chunk of subruns in $2...\n"

    runNumber=$2
    lowSubrun=$3
    highSubrun=$4

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

    ./bin/sumChunk "$analysisDirectoryName" "$experiment" "$runNumber" "$lowSubrun" "$highSubrun"

fi
