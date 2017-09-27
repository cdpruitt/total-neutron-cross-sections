#ifndef PLOTTING_CONSTANTS_H
#define PLOTTING_CONSTANTS_H

#include <math.h>

#include "physicalConstants.h"
#include "runSpecificConstants.h"
#include <vector>
#include <string>

/*****************************************************************************/
/* Plotting constants*/

const int TOF_LOWER_BOUND = 0;   // lower bound of cross-section plots
const int TOF_UPPER_BOUND = ceil(MICRO_LENGTH); // upper bound of cross-section plots
const int TOF_RANGE = TOF_UPPER_BOUND-TOF_LOWER_BOUND;
const int TOF_BINS = TOF_RANGE*4; // for plots with time units as abscissa

const double ENERGY_LOWER_BOUND = 1.7;   // lower bound of cross-section plots, in MeV
const double ENERGY_UPPER_BOUND = 800; // upper bound of cross-section plots, in MeV
const unsigned int NUMBER_ENERGY_BINS = 300;    // for plots with energy units as abscissa

// for RKEToTOF functions

const unsigned int NUMBER_TOF_BINS = 300;    // for plots with energy units as abscissa

const unsigned int ENERGY_RANGE = 300; // upper bound of time-of-flight plots, in ns
const unsigned int ENERGY_BINS = 300; // for plots with time units as abscissa


#endif /* PLOTTING_CONSTANT_H */
