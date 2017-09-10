#ifndef ANALYSIS_CONSTANTS_H
#define ANALYSIS_CONSTANTS_H

#include <vector>
#include <string>
#include <utility>

/*****************************************************************************/
/* Constants used to structure and analyze experimental data */

/******************************************************************************/
// "Sync window" defines the acceptable time variation between consecutive macropulses
// (with respect to the macropulse period) in order for those macropulses to be
// considered still "in sync" with the macropulse period
const double SYNC_WINDOW = 0.04; // as fraction of MACRO_PERIOD 

// "Target changer charge gates" are used to assign the target changer position
// based on the target changer signal's integrated charge
const std::vector<std::pair<int,int>> tarGates = {
    std::pair<int,int> (50,2500), // position 0 gates
    std::pair<int,int> (5000,10000), // position 1 gates
    std::pair<int,int> (12000,17000), // position 2 gates
    std::pair<int,int> (18000,23000), // position 3 gates
    std::pair<int,int> (25000,30000), // position 4 gates
    std::pair<int,int> (31000,36000), // position 5 gates
    std::pair<int,int> (37000,43000) // position 6 gates
};

const double SCALEDOWN = 10; // (for debugging) only sort (total/SCALEDOWN) events

// experimentally-determined digitizer deadtime
const int DEADTIME_PERIOD = 225; // in ns
const int DEADTIME_TRANSITION_PERIOD = 8; // in ns

const unsigned int TARGET_CHANGER_LED_THRESHOLD = 1500; // in ADC units
const unsigned int MAIN_DETECTOR_LED_THRESHOLD = 3000; // in ADC units

const unsigned int TIME_CHECK_TOLERANCE = 10; // in ns

const std::vector<int> tcFineTimeThresholds = {1000, 1750, 2650, 3600, 4400, 5250}; // in ADC units

// Indicate the range of times considered to be gamma rays (for the purposes of
// counting gamma rays)
const double GAMMA_WINDOW[2] = {82,88};

const std::vector<std::string> positionNames = {"blank", "target1", "target2", "target3", "target4", "target5"};

const std::vector<std::string> targetNamesWaveform = {"blankWaveform", "shortCarbonWaveform", "longCarbonWaveform", "Sn112Waveform", "NatSnWaveform", "Sn124Waveform"}; 

const std::vector<std::string> detectorNames = {"summedDet"};
//const std::vector<std::string> detectorNames = {"highThresholdDet","lowThresholdDet"};

const std::vector<double> manualTimeOffsets = {0, -0.7, -1.04, -1.23, -1.39, -1.44};

#endif
