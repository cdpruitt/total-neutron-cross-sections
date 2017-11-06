#ifndef FILL_ADVANCED_HISTOS_H
#define FILL_ADVANCED_HISTOS_H

#include "TH1D.h"
#include "plots.h"

#include "../include/GammaCorrection.h"

int fillCSHistos(std::string inputFileName, std::ofstream& log, std::string treeName, std::vector<GammaCorrection> gammaCorrectionList, std::string outputFileName);

TH1D* convertTOFtoEnergy(TH1D* tof, std::string name);

#endif /* FILL_ADVANCED_HISTOS_H */
