#ifndef FILL_ADVANCED_HISTOS_H
#define FILL_ADVANCED_HISTOS_H

#include "TH1D.h"
#include "plots.h"

int fillAdvancedHistos(std::string inputFileName, std::string treeName, std::string outputFileName);

TH1D* convertTOFtoEnergy(TH1D* tof, std::string name);

#endif /* FILL_ADVANCED_HISTOS_H */
