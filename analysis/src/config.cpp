#include "../include/config.h"
#include "../include/experiment.h"

#include <string>
#include <vector>
#include <utility> // for "pair"

using namespace std;

FacilityConfig::FacilityConfig(std::vector<std::string> facilityConfig)
{
    if(facilityConfig.size()!=6)
    {
        std::cerr << "Error: tried to create facility configuration for current run, but " << facilityConfig.size() << " facility parameters were read in instead of the correct value of 6." << std::endl;
        exit(1);
    }

    MACRO_FREQUENCY = stod(facilityConfig[0]);
    MICROS_PER_MACRO = stoi(facilityConfig[1]);
    MICRO_LENGTH = stod(facilityConfig[2]);

    FLIGHT_DISTANCE = stod(facilityConfig[3]);

    FIRST_GOOD_MICRO = stod(facilityConfig[4]);
    LAST_GOOD_MICRO = stod(facilityConfig[5]);
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
}

TimeConfig::TimeConfig(std::vector<std::string> timeConfig)
{
    if(timeConfig.size()!=9)
    {
        std::cerr << "Error: tried to create time offsets configuration for current run, but " << timeConfig.size() << " time offsets were read in" << std::endl;
        exit(1);
    }

    for(int i=0; i<timeConfig.size()-1; i++)
    {
        offsets.push_back(stod(timeConfig[i]));
    }

    GAMMA_WINDOW_SIZE = stod(timeConfig[8]);
}

DigitizerConfig::DigitizerConfig(vector<pair<int, string>> channelMap, vector<string> digitizerConfig)
{
    if(channelMap.size()==0)
    {
        cerr << "Error: tried to create digitizer configuration for current run, but no channel data could be read in." << endl;
        exit(1);
    }

    CHANNEL_MAP = channelMap;

    if(digitizerConfig.size()==0)
    {
        cerr << "Error: tried to create digitizer configuration for current run, but no config data were read in." << endl;
        exit(1);
    }

    SAMPLE_PERIOD = stoi(digitizerConfig[0]);
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
}

CSConfig::CSConfig(std::vector<std::string> v)
{
    if(v.size() >= 8)
    {
        std::cerr << "Error: tried to create plot configuration for current run, but wrong number of parameters were read in." << std::endl;
        exit(1);
    }

    DETECTOR_NAMES = v;
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
}

Config::Config(std::string expName, int runNumber)
{
    facility = readFacilityConfig(expName, runNumber);
    softwareCFD = readSoftwareCFDConfig(expName, runNumber);
    target = readTargetConfig(expName, runNumber);
    cs = readCSConfig(expName, runNumber);
    plot = readPlotConfig(expName, runNumber);
    digitizer = readDigitizerConfig(expName, runNumber);
    time = readTimeConfig(expName, runNumber);
    deadtime = readDeadtimeConfig(expName);
    analysis = readAnalysisConfig(expName);

    std::cout << "Finished reading experiment config data." << std::endl;
}
