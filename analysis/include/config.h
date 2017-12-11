#ifndef EXPERIMENTAL_CONFIG_H
#define EXPERIMENTAL_CONFIG_H

#include <string>
#include <vector>
#include <iostream>

struct AnalysisConfig
{
    public:
        AnalysisConfig() {}

        std::string RAW_TREE_FILE_NAME;
        std::string MACROPULSE_ASSIGNED_FILE_NAME;
        std::string PASSED_VETO_FILE_NAME;
        std::string HISTOGRAM_FILE_NAME;
        std::string ENERGY_PLOTS_FILE_NAME;

        std::string DPP_TREE_NAME;
        std::string WAVEFORM_TREE_NAME;
        std::string MACROPULSE_TREE_NAME;
        std::string GAMMA_CORRECTION_TREE_NAME;
};

struct DeadtimeConfig
{
    public:
        DeadtimeConfig() {}

        double LOGISTIC_K;
        double LOGISTIC_MU;
};

struct FacilityConfig
{
    public:
        FacilityConfig() {}
        FacilityConfig(std::vector<std::string> facilityConfig);

        double MACRO_FREQUENCY;   // frequency of macropulses, in Hz
        unsigned int MICROS_PER_MACRO;  // number of micropulses in each macropulse
        double MICRO_LENGTH;      // micropulse duration, in ns
        double FLIGHT_DISTANCE;   // detector distance from neutron source, in cm

        unsigned int FIRST_GOOD_MICRO; // first micropulse to be used for cs calculation
        unsigned int LAST_GOOD_MICRO;  // last micropulse to be used for cs calculation
};

struct SoftwareCFDConfig
{
    public:
        SoftwareCFDConfig() {}
        SoftwareCFDConfig(std::vector<std::string> softwareCFDConfig);

        double CFD_FRACTION;    // allowed range: 0 to 1
        double CFD_DELAY; // in samples
        double CFD_ZC_TRIGGER_THRESHOLD;
        // prevent ZC triggers until the sum of the delayed, original waveform and
        // the reduced-amplitude, opposite-sign waveform sum crosses a certain
        // threshold. This largely excludes spurious ZC crossings from noise before
        // the main pulse of the original waveform

        double CFD_TIME_OFFSET; // offsets the calculated CFD time so that the average CFD time correction is 0 (in samples)
};

struct TimeConfig
{
    public:
        TimeConfig() {}
        TimeConfig(std::vector<std::string> timeConfig);

        std::vector<double> offsets;

        double GAMMA_WINDOW_SIZE; // width of window used to identify gammas for gamma time correction to macropulse times
};

struct DigitizerConfig
{
    public:
        DigitizerConfig() {}
        DigitizerConfig(std::vector<std::pair<unsigned int, std::string>> channelMap, std::vector<std::string> digitizerConfig);

        std::vector<std::pair<unsigned int, std::string>> CHANNEL_MAP;
        double SAMPLE_PERIOD;
};

struct TargetConfig
{
    public:
        TargetConfig() {}
        TargetConfig(std::vector<std::string> targetPositions, std::vector<std::pair<int,int>> targetChangerGates);

        std::vector<std::pair<int,int>> TARGET_GATES;
        std::vector<std::string> TARGET_ORDER;
};

struct CSConfig
{
    public:
        CSConfig() {}
        CSConfig(std::vector<std::string> v);
    
        std::vector<std::string> DETECTOR_NAMES;
};

struct PlotConfig
{
    public:
        PlotConfig() {}
        PlotConfig(std::vector<std::string> v);

        unsigned int TOF_LOWER_BOUND;   // lower bound of cross-section plots
        unsigned int TOF_UPPER_BOUND; // upper bound of cross-section plots
        unsigned int TOF_BINS_PER_NS;
        unsigned int TOF_RANGE;
        unsigned int TOF_BINS;

        double ENERGY_LOWER_BOUND;   // lower bound of cross-section plots, in MeV
        double ENERGY_UPPER_BOUND; // upper bound of cross-section plots, in MeV
        unsigned int NUMBER_ENERGY_BINS;    // for plots with energy units as abscissa
};

struct Config
{
    public:
        Config() {}
        Config(std::string expName, int runNumber);

        AnalysisConfig analysis;
        DeadtimeConfig deadtime;
        FacilityConfig facility;
        SoftwareCFDConfig softwareCFD;
        TargetConfig target;
        CSConfig cs;
        PlotConfig plot;
        DigitizerConfig digitizer;
        TimeConfig time;
};

extern Config config;

#endif /* EXPERIMENTAL_CONFIG_H */
