#ifndef PLOTTING_CONSTANTS_H
#define PLOTTING_CONSTANTS_H

#include "physicalConstants.h"
#include "analysisConstants.h"
#include <vector>
#include <string>

/*****************************************************************************/
/* Plotting constants*/

const double ENERGY_LOWER_BOUND = 1;   // lower bound of cross-section plots, in MeV
const double ENERGY_UPPER_BOUND = 600; // upper bound of cross-section plots, in MeV
const int NUMBER_ENERGY_BINS = 30;    // for plots with energy units as abscissa

const int TOF_RANGE = MICRO_LENGTH+1; // upper bound of time-of-flight plots, in ns
const int TOF_BINS = 3600; // for plots with time units as abscissa

#endif
