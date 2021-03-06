#ifndef FILL_ADVANCED_HISTOS_H
#define FILL_ADVANCED_HISTOS_H

#include "TH1D.h"
#include "plots.h"

#include "../include/GammaCorrection.h"

int fillCSHistos(std::string vetoedInputFileName, std::string nonVetoInputFileName, bool useVetoPaddle, std::string macropulseFileName, std::string gammaCorrectionFileName, std::ofstream& log, std::string outputFileName);

int fillMonitorHistos(std::string inputFileName, std::string macropulseFileName, std::ofstream& log, std::string outputFileName);

TH1D* convertTOFtoEnergy(TH1D* tof, std::string name);

#endif /* FILL_ADVANCED_HISTOS_H */
