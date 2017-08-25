#ifndef HISTOS_H
#define HISTOS_H

#include "TH1D.h"
#include "plots.h"

int histos(std::string sortedFileName, std::string vetoedFileName, std::string histoFileName, std::vector<std::string> channelMap);
TH1D* convertTOFtoEnergy(TH1D* tof, std::string name);

#endif
