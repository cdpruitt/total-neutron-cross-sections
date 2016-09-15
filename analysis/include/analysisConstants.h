#ifndef ANALYSIS_CONSTANTS_H
#define ANALYSIS_CONSTANTS_H

#include <vector>
#include <string>

/*****************************************************************************/
/* Constants used to structure and analyze experimental data */

/******************************************************************************/
// "Sync window" defines the acceptable time variation between consecutive macropulses
// (with respect to the macropulse period) in order for those macropulses to be
// considered still "in sync" with the macropulse period
const double SYNC_WINDOW = 0.005; // as fraction of MACRO_PERIOD 

// "Target changer charge gates" are used to assign the target changer position
// based on the target changer signal's integrated charge
const int tarGate[12] = {5000,7500,11000,13500,17000,20000,23000,27000,29000,33000,35000,39000};
            // Position: 1low  1hi  2low   2hi  3low   3hi  4low   4hi  5low   5hi  6low   6hi

// Establish which channels are active in each mode
const std::vector<std::string> activeDPPChannels = {"ch0","ch2","ch4","ch6"};
const std::vector<std::string> activeWaveformChannels = {"ch0","ch2","ch4"};

const double SCALEDOWN = 1; // (for debugging) only sort (total/SCALEDOWN) events

// experimentally-determined  digitizer deadtime
const int DEADTIME_PERIOD = 183;

// Indicate the range of times considered to be gamma rays (for the purposes of
// counting gamma rays)
const double GAMMA_WINDOW[2] = {80,100};

const std::vector<std::string> targetNamesWaveform = {"blankWaveform", "shortCarbonWaveform", "longCarbonWaveform", "Sn112Waveform", "NatSnWaveform", "Sn124Waveform"}; 


#endif
