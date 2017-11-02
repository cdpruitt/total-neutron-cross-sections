#ifndef EXPERIMENTAL_CONFIG_H
#define EXPERIMENTAL_CONFIG_H

#include <string>
#include <vector>
#include <iostream>

struct FacilityConfig
{
    public:
        FacilityConfig() {}
        FacilityConfig(std::vector<std::string> facilityConfig);

        double MACRO_FREQUENCY;   // frequency of macropulses, in Hz
        unsigned int MICROS_PER_MACRO;  // number of micropulses in each macropulse
        double MICRO_LENGTH;      // micropulse duration, in ns
        double FLIGHT_DISTANCE;   // detector distance from neutron source, in cm
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

struct TimeOffsetsConfig
{
    public:
        TimeOffsetsConfig() {}
        TimeOffsetsConfig(std::vector<std::string> timeOffsetsConfig);

        double TARGET_CHANGER_TIME_OFFSET; // timing delay between macropulse start and target changer position channel (in ns)

        double MONITOR_TIME_OFFSET; // timing delay between macropulse start and monitor channel (in ns)

        double DETECTOR_TIME_OFFSET; // timing delay between macropulse start and main detector channel (in ns)

        double VETO_TIME_OFFSET; // timing delay between macropulse start and veto detector channel (in ns)

        double HIGH_T_DET_TIME_OFFSET; // timing delay between macropulse start and high threshold detector channel (in ns)

        double GAMMA_WINDOW_SIZE; // width of window used to identify gammas for gamma time correction to macropulse times
};

struct DigitizerConfig
{
    DigitizerConfig() {}
    double SAMPLE_PERIOD = 2; // digitizer sample rate, in ns
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

        FacilityConfig facility;
        SoftwareCFDConfig softwareCFD;
        TargetConfig target;
        CSConfig cs;
        PlotConfig plot;
        DigitizerConfig digitizer;
        TimeOffsetsConfig timeOffsets;
};

extern Config config;

#endif /* EXPERIMENTAL_CONFIG_H */
