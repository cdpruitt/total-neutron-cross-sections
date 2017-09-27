#ifndef EXPERIMENTAL_CONFIG_H
#define EXPERIMENTAL_CONFIG_H

struct FacilityConfig
{
    double MACRO_FREQUENCY = 120;   // frequency of beam macropulses, in Hz
    double MACRO_LENGTH = 625000;   // macropulse duration, in ns
    double MICRO_LENGTH = 1788.814; // micropulse duration, in ns
    double FLIGHT_DISTANCE = 2710;  // detector distance from neutron
    // source, in cm
}

struct TimeConfig
{
    double CFD_FRACTION = 0.75;
    unsigned int CFD_DELAY = 2;
    double CFD_ZC_TRIGGER_THRESHOLD = 30*CFD_FRACTION;

    double TC_CFD_FRACTION = 0.5;
    unsigned int TC_CFD_DELAY = 2;
    double TC_ZC_TRIGGER_THRESHOLD = 500*TC_CFD_FRACTION;

    double TC_FINETIME_THRESHOLD = 3000;

    const double MACROPULSE_OFFSET = 775.35; // timing delay between the digitizer (channel 0)
    // and the facility's RF clock (due to cable delay,
    // NIM logic, etc.), in ns
    const double VETO_OFFSET = 8;         // timing delay of the veto paddle after the main
    // detector channel, in ns

    const std::vector<int> tcFineTimeThresholds = {1000, 1750, 2650, 3600, 4400, 5250}; // in ADC units

    const std::vector<double> MACROTIME_TARGET_DRIFT = {0, 0.30, 0.42, 0.52, 0.58, 0.66};

    const unsigned int TARGET_CHANGER_LED_THRESHOLD = 1500; // in ADC units
    const unsigned int MAIN_DETECTOR_LED_THRESHOLD = 3000; // in ADC units
}

struct TargetConfig
{
}

struct CSConfig
{
    const std::vector<std::string> detectorNames = {"summedDet","highTDet"};
}

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
}

#endif /* EXPERIMENTAL_CONFIG_H */
