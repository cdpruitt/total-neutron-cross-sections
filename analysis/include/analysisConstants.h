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
const double SYNC_WINDOW = 0.02; // as fraction of MACRO_PERIOD 

// "Target changer charge gates" are used to assign the target changer position
// based on the target changer signal's integrated charge
const std::vector<std::pair<int,int>> tarGates = {
    std::pair<int,int> (1500,2500), // position 1 gates
    std::pair<int,int> (3000,4400), // position 2 gates
    std::pair<int,int> (4600,6000), // position 3 gates
    std::pair<int,int> (6500,7800), // position 4 gates
    std::pair<int,int> (8000,9500), // position 5 gates
    std::pair<int,int> (9600,11500) // position 6 gates
};

const double SCALEDOWN = 1; // (for debugging) only sort (total/SCALEDOWN) events

// experimentally-determined  digitizer deadtime
const int DEADTIME_PERIOD = 152;

// Indicate the range of times considered to be gamma rays (for the purposes of
// counting gamma rays)
const double GAMMA_WINDOW[2] = {75,95};

const std::vector<std::string> positionNames = {"blank", "target1", "target2", "target3", "target4", "target5"};

const std::vector<std::string> targetNamesWaveform = {"blankWaveform", "shortCarbonWaveform", "longCarbonWaveform", "Sn112Waveform", "NatSnWaveform", "Sn124Waveform"}; 

const std::vector<std::string> detectorNames = {"highThresholdDet","lowThresholdDet"};

#endif
