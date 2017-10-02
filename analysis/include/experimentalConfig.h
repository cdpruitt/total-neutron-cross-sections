#ifndef EXPERIMENTAL_CONFIG_H
#define EXPERIMENTAL_CONFIG_H

struct FacilityConfig
{
    double MACRO_FREQUENCY = 120;   // frequency of beam macropulses, in Hz
    double MACRO_LENGTH = 625000;   // macropulse duration, in ns
    double MICRO_LENGTH = 1788.814; // micropulse duration, in ns
    double FLIGHT_DISTANCE = 2710;  // detector distance from neutron
    // source, in cm
};

struct TimeConfig
{
    double CFD_FRACTION = 0.75;
    unsigned int CFD_DELAY = 3;
    double CFD_ZC_TRIGGER_THRESHOLD = 50*CFD_FRACTION;

    double TC_CFD_FRACTION = 0.5;
    unsigned int TC_CFD_DELAY = 2;
    double TC_ZC_TRIGGER_THRESHOLD = 500*TC_CFD_FRACTION;

    double FINE_TIME_OFFSET = 20.96; // in samples

    double MACROPULSE_FINE_TIME_THRESHOLD = 11500;

    double MACROPULSE_TARGET_TIME_DIFFERENCE = 68; // in ns;

    const double SUMMED_DETECTOR_TIME_OFFSET = 775.35; // timing delay between the digitizer (channel 0)
    // and the facility's RF clock (due to cable delay,
    // NIM logic, etc.), in ns

    const double MONITOR_DETECTOR_TIME_OFFSET = 775.35; // timing delay between the digitizer (channel 0)
    // and the facility's RF clock (due to cable delay,
    // NIM logic, etc.), in ns

    const double VETO_OFFSET = 8;         // timing delay of the veto paddle after the main
    // detector channel, in ns

    const std::vector<int> tcFineTimeThresholds = {1000, 1750, 2650, 3600, 4400, 5250}; // in ADC units

    const std::vector<double> MACROTIME_TARGET_DRIFT = {0, 0.30, 0.42, 0.52, 0.58, 0.66};

    const unsigned int TARGET_CHANGER_LED_THRESHOLD = 1500; // in ADC units
    const unsigned int MAIN_DETECTOR_LED_THRESHOLD = 3000; // in ADC units

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
    const std::vector<std::string> DETECTOR_NAMES = {"summedDet"};
};

struct PlotConfig
{
    const int TOF_LOWER_BOUND = 0;   // lower bound of cross-section plots
    const int TOF_UPPER_BOUND = 1789; // upper bound of cross-section plots
    const int TOF_RANGE = TOF_UPPER_BOUND-TOF_LOWER_BOUND;
    const int TOF_BINS = TOF_RANGE*4; // for plots with time units as abscissa

    const double ENERGY_LOWER_BOUND = 1.7;   // lower bound of cross-section plots, in MeV
    const double ENERGY_UPPER_BOUND = 800; // upper bound of cross-section plots, in MeV
    const unsigned int NUMBER_ENERGY_BINS = 300;    // for plots with energy units as abscissa

    // for RKEToTOF functions

    const unsigned int NUMBER_TOF_BINS = 300;    // for plots with energy units as abscissa

    const unsigned int ENERGY_RANGE = 300; // upper bound of time-of-flight plots, in ns
    const unsigned int ENERGY_BINS = 300; // for plots with time units as abscissa
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
