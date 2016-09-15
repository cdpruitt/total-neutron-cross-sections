#ifndef PLOTTING_CONSTANTS_H
#define PLOTTING_CONSTANTS_H

#include <vector>
#include <string>

/*****************************************************************************/
/* Plotting constants*/

const double ENERGY_LOWER_BOUND = 4;   // lower bound of cross-section plots, in MeV
const double ENERGY_UPPER_BOUND = 600; // upper bound of cross-section plots, in MeV
const int NUMBER_ENERGY_BINS = 200;    // for plots with energy units as abscissa

const int TOF_RANGE = 1800; // upper bound of time-of-flight plots, in ns
const int TOF_BINS = 18000; // for plots with time units as abscissa

const std::vector<std::string> positionNames =
{"blank", "target1", "target2", "target3", "target4", "target5"}; 
// segregate histograms by target position, with labels given by 'positionNames'

const std::string dirs[4] = {"targetChanger","monitor","detS","scavenger"};
// segregate directories by data channel, with labels given by 'dirs'

#endif
