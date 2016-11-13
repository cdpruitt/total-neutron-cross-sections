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
const double SYNC_WINDOW = 0.02; // as fraction of MACRO_PERIOD 

// "Target changer charge gates" are used to assign the target changer position
// based on the target changer signal's integrated charge
const int tarGate[12] = {1500,2500,3000,4400,4600,6000,6500,7800,8000,9500,9600,11500};
            // Position: 1low  1hi 2low  2hi 3low  3hi 4low 4hi  5low 5hi  6low   6hi

// Establish which channels are active, and map them to function
const std::vector<std::pair<std::string,std::string>> channelMap = {
                  std::pair<std::string,std::string>("ch0","targetChanger"),
                  std::pair<std::string,std::string>("ch1","ch1"),
                  std::pair<std::string,std::string>("ch2","monitor"),
                  std::pair<std::string,std::string>("ch3","ch3"),
                  std::pair<std::string,std::string>("ch4","highThresholdDet"),
                  std::pair<std::string,std::string>("ch5","lowThresholdDet"),
                  std::pair<std::string,std::string>("ch6","vetoPaddle"),
                  std::pair<std::string,std::string>("ch7","ch7"),
};

const double SCALEDOWN = 1; // (for debugging) only sort (total/SCALEDOWN) events

// experimentally-determined  digitizer deadtime
const int DEADTIME_PERIOD = 152;

// Indicate the range of times considered to be gamma rays (for the purposes of
// counting gamma rays)
const double GAMMA_WINDOW[2] = {75,95};

const std::vector<std::string> targetNamesWaveform = {"blankWaveform", "shortCarbonWaveform", "longCarbonWaveform", "Sn112Waveform", "NatSnWaveform", "Sn124Waveform"}; 


#endif
