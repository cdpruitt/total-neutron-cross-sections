#ifndef FILL_BASIC_HISTOS_H
#define FILL_BASIC_HISTOS_H

#include "TH1D.h"
#include "plots.h"

#include "../include/GammaCorrection.h"

int fillBasicHistos(std::string inputFileName, std::ofstream& log, std::string treeName, std::vector<GammaCorrection> gammaCorrectionList, std::string outputFileName);

#endif /* FILL_BASIC_HISTOS_H */
