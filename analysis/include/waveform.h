#ifndef WAVEFORM_H
#define WAVEFORM_H

#include <vector>

#include "plots.h"

void waveform(std::string inFileName, TFile*& outFile, const std::vector<Target*>& targets, std::vector<Plots*>& plots);

#endif
