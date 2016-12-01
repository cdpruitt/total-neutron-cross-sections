#ifndef HISTOS_H
#define HISTOS_H

#include "plots.h"

int histos(std::string sortedFileName, std::string vetoedFileName, std::string histoFileName, std::vector<std::string> channelMap);
void correctForDeadtime(std::string histoFileName, std::string waveformFileName);

#endif
