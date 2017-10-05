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

read -r experiment<experiment.txt # find out which experimental dataset to analyze

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

./bin/sumAll "$analysisDirectoryName" "$experiment"
