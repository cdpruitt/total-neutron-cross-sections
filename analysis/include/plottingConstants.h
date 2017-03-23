#ifndef PLOTTING_CONSTANTS_H
#define PLOTTING_CONSTANTS_H

#include "physicalConstants.h"
#include "analysisConstants.h"
#include <vector>
#include <string>

/*****************************************************************************/
/* Plotting constants*/

const double TOF_LOWER_BOUND = 0;   // lower bound of cross-section plots
const double TOF_UPPER_BOUND = MICRO_LENGTH; // upper bound of cross-section plots
const int TOF_RANGE = TOF_UPPER_BOUND-TOF_LOWER_BOUND;
const int TOF_BINS = 3600; // for plots with time units as abscissa

const double ENERGY_LOWER_BOUND = 1.7;   // lower bound of cross-section plots, in MeV
const double ENERGY_UPPER_BOUND = 800; // upper bound of cross-section plots, in MeV
const int NUMBER_ENERGY_BINS = 300;    // for plots with energy units as abscissa

// for RKEToTOF functions

const int NUMBER_TOF_BINS = 300;    // for plots with energy units as abscissa

const int ENERGY_RANGE = 300; // upper bound of time-of-flight plots, in ns
const int ENERGY_BINS = 300; // for plots with time units as abscissa


#endif
