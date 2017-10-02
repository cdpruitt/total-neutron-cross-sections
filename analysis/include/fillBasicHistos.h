#ifndef FILL_BASIC_HISTOS_H
#define FILL_BASIC_HISTOS_H

#include "TH1D.h"
#include "plots.h"

int fillBasicHistos(std::string inputFileName, std::string treeName, std::string outputFileName);
TH1D* convertTOFtoEnergy(TH1D* tof, std::string name);

#endif /* FILL_BASIC_HISTOS_H */
