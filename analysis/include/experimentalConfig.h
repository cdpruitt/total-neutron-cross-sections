#ifndef EXPERIMENTAL_CONFIG_H
#define EXPERIMENTAL_CONFIG_H

struct FacilityConfig
{
    double MACRO_FREQUENCY;   // frequency of macropulses, in Hz
    double MACRO_LENGTH;      // macropulse duration, in ns
    double MICRO_LENGTH;      // micropulse duration, in ns
    double FLIGHT_DISTANCE;   // detector distance from neutron source, in cm
};

struct SoftwareCFDConfig
{
    double CFD_FRACTION;    // allowed range: 0 to 1

    unsigned int CFD_DELAY; // in samples

    double CFD_ZC_TRIGGER_THRESHOLD;
    // prevent ZC triggers until the sum of the delayed, original waveform and
    // the reduced-amplitude, opposite-sign waveform sum crosses a certain
    // threshold. This largely excludes spurious ZC crossings from noise before
    // the main pulse of the original waveform

    double CFD_TIME_OFFSET; // offsets the calculated CFD time so that the average CFD time correction is 0 (in samples)
}

struct TimeOffsetsConfig
{
    double TARGET_CHANGER_TIME_OFFSET; // timing delay between macropulse start and target changer position channel (in ns)

    double MONITOR_TIME_OFFSET; // timing delay between macropulse start and monitor channel (in ns)

    double DETECTOR_TIME_OFFSET; // timing delay between macropulse start and main detector channel (in ns)

    double VETO_TIME_OFFSET; // timing delay between macropulse start and veto detector channel (in ns)
}

struct DigitizerConfig
{
    const double SAMPLE_PERIOD = 2; // digitizer sample rate, in ns
};

struct TargetConfig
{
    // "Target changer charge gates" are used to assign the target changer position
    // based on the target changer signal's integrated charge
    const std::vector<std::pair<int,int>> TARGET_GATES = {
        std::pair<int,int> (50,2500), // position 0 gates
        std::pair<int,int> (5000,10000), // position 1 gates
        std::pair<int,int> (12000,17000), // position 2 gates
        std::pair<int,int> (18000,23000), // position 3 gates
        std::pair<int,int> (25000,30000), // position 4 gates
        std::pair<int,int> (31000,36000), // position 5 gates
        std::pair<int,int> (37000,43000) // position 6 gates
    };

    const std::vector<std::string> POSITION_NAMES = {"target0", "target1", "target2", "target3", "target4", "target5"};
};

struct CSConfig
{
    std::vector<std::string> DETECTOR_NAMES;
};

struct PlotConfig
{
    unsigned int TOF_LOWER_BOUND;   // lower bound of cross-section plots
    unsigned int TOF_UPPER_BOUND; // upper bound of cross-section plots
    unsigned int TOF_RANGE = TOF_UPPER_BOUND-TOF_LOWER_BOUND;
    unsigned int TOF_BINS = TOF_RANGE*4; // for plots with time units as abscissa

    double ENERGY_LOWER_BOUND;   // lower bound of cross-section plots, in MeV
    double ENERGY_UPPER_BOUND; // upper bound of cross-section plots, in MeV
    unsigned int NUMBER_ENERGY_BINS;    // for plots with energy units as abscissa
};

struct ExperimentalConfig
{
    // facility parameters
    FacilityConfig facilityConfig;

    // time calculation parameters
    TimeConfig timeConfig;

    // target parameters
    TargetConfig targetConfig;

    // cross section calculation parameters
    CSConfig csConfig;

    // plotting parameters
    PlotConfig plotConfig;

};

extern ExperimentalConfig experimentalConfig;

#endif /* EXPERIMENTAL_CONFIG_H */
