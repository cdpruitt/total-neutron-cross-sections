#ifndef CORRECT_FOR_DEADTIME_H
#define CORRECT_FOR_DEADTIME_H

#include "GammaCorrection.h"

#include "TH1.h"

#include <fstream>
#include <string>
#include <vector>

int generateDeadtimeCorrection(std::string inputFileName, std::ofstream& log, std::string outputFileName);

int applyDeadtimeCorrection(std::string inputFileName, std::string deadtimeHistoFileName, std::string macroFileName, std::string gammaCorrectionFileName, std::ofstream& log, std::string outputFileName);

#endif /* CORRECT_FOR_DEADTIME_H */
