#ifndef CORRECT_FOR_DEADTIME_H
#define CORRECT_FOR_DEADTIME_H

#include <vector>

#include "TH1.h"

int generateDeadtimeCorrection(TH1D* tof,
        unsigned int numberOfMacros, std::vector<double>& deadtimeCorrectionList);

void applyDeadtimeCorrection(TH1D* rawTOF, TH1D* correctedTOF, const std::vector<double>& deadtimesPerBin);

#endif /* CORRECT_FOR_DEADTIME_H */
