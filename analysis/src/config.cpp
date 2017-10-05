#include "../include/config.h"
#include "../include/experiment.h"

#include <string>
#include <vector>
#include <utility> // for "pair"

using namespace std;

FacilityConfig::FacilityConfig(std::vector<std::string> facilityConfig)
{
    if(facilityConfig.size()!=4)
    {
        std::cerr << "Error: tried to create facility configuration for current run, but " << facilityConfig.size() << " facility parameters were read in instead of the correct value of 4." << std::endl;
        exit(1);
    }

    MACRO_FREQUENCY = stod(facilityConfig[0]);
    MACRO_LENGTH = stod(facilityConfig[1]);
    MICRO_LENGTH = stod(facilityConfig[2]);
    FLIGHT_DISTANCE = stod(facilityConfig[3]);

    std::cout << "Successfully read facility config data." << std::endl;
}

SoftwareCFDConfig::SoftwareCFDConfig(std::vector<std::string> softwareCFDConfig)
{
    if(softwareCFDConfig.size()!=4)
    {
        std::cerr << "Error: tried to create softwareCFD configuration for current run, but " << softwareCFDConfig.size() << " software CFD parameters were read in instead of the correct value of 4." << std::endl;
        exit(1);
    }

    CFD_FRACTION = stod(softwareCFDConfig[0]);
    CFD_DELAY = stod(softwareCFDConfig[1]);
    CFD_ZC_TRIGGER_THRESHOLD = stod(softwareCFDConfig[2]);
    CFD_TIME_OFFSET = stod(softwareCFDConfig[3]);

    std::cout << "Successfully read software CFD config data." << std::endl;
}

TimeOffsetsConfig::TimeOffsetsConfig(std::vector<std::string> timeOffsetsConfig)
{
    if(timeOffsetsConfig.size()!=4)
    {
        std::cerr << "Error: tried to create time offsets configuration for current run, but " << timeOffsetsConfig.size() << " time offsets were read in instead of the correct value of 4." << std::endl;
        exit(1);
    }

    TARGET_CHANGER_TIME_OFFSET = stod(timeOffsetsConfig[0]);
    MONITOR_TIME_OFFSET = stod(timeOffsetsConfig[1]);
    DETECTOR_TIME_OFFSET = stod(timeOffsetsConfig[2]);
    VETO_TIME_OFFSET = stod(timeOffsetsConfig[3]);

    std::cout << "Successfully read time offset config data." << std::endl;
}

TargetConfig::TargetConfig(std::vector<std::string> targetPositions, std::vector<std::pair<int,int>> targetChangerGates)
{
    if(targetPositions.size()==0)
    {
        std::cerr << "Error: tried to create target configuration for current run, but zero target positions were read in." << std::endl;
        exit(1);
    }

    TARGET_ORDER = targetPositions;

    if(targetChangerGates.size()==0)
    {
        std::cerr << "Error: tried to create target configuration for current run, but no target changer gates were read in." << std::endl;
        exit(1);
    }

    TARGET_GATES = targetChangerGates;

    std::cout << "Successfully read target changer config data." << std::endl;
}

CSConfig::CSConfig(std::vector<std::string> v)
{
    if(v.size() != 1)
    {
        std::cerr << "Error: tried to create plot configuration for current run, but wrong number of parameters were read in." << std::endl;
        exit(1);
    }

    DETECTOR_NAMES = std::vector<std::string>();

    for(std::string s : v)
    {
        DETECTOR_NAMES.push_back(s);
    }

    std::cout << "Successfully read cross section config data." << std::endl;
}

PlotConfig::PlotConfig(std::vector<std::string> v)
{
    if(v.size() != 6)
    {
        std::cerr << "Error: tried to create plot configuration for current run, but wrong number of parameters were read in." << std::endl;
        exit(1);
    }

    TOF_LOWER_BOUND = stoi(v[0]);
    TOF_UPPER_BOUND = stoi(v[1]);
    TOF_BINS_PER_NS = stoi(v[2]);
    TOF_RANGE = TOF_UPPER_BOUND-TOF_LOWER_BOUND;
    TOF_BINS = TOF_RANGE*TOF_BINS_PER_NS;

    ENERGY_LOWER_BOUND = stod(v[3]);
    ENERGY_UPPER_BOUND = stod(v[4]);
    NUMBER_ENERGY_BINS = stoi(v[5]);

    std::cout << "Successfully read plot config data." << std::endl;
}

Config::Config(std::string expName, int runNumber)
{
    facilityConfig = readFacilityConfig(expName, runNumber);
    softwareCFDConfig = readSoftwareCFDConfig(expName, runNumber);
    targetConfig = readTargetConfig(expName, runNumber);
    csConfig = readCSConfig(expName, runNumber);
    plotConfig = readPlotConfig(expName, runNumber);
    digitizerConfig = DigitizerConfig();
    timeOffsetsConfig = readTimeOffsetsConfig(expName, runNumber);

    std::cout << "Finished reading experiment config data." << std::endl;
}

