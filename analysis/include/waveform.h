#ifndef WAVEFORM_H
#define WAVEFORM_H

#include <vector>

#include "plots.h"

void waveform(std::string inFileName,std::string outFileName,std::string experiment);

void calculateDeadtime(std::vector<long> microsPerTarget, std::vector<Plots*>& plots);
#endif
