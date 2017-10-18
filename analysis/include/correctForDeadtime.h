#ifndef CORRECT_FOR_DEADTIME_H
#define CORRECT_FOR_DEADTIME_H

#include <vector>

#include "TH1.h"
#include "TDirectory.h"

int correctForDeadtime(TH1D*& rawTOF, TH1D*& correctedTOF, const double deadtimeBins, const double deadtimeTransitionBins, const unsigned int& numberOfMacros);

int correctForDeadtime2(TH1D*& TOFtoCorrect, TH1D*& correctedTOF, const double deadtimeBins, const double deadtimeTransitionBins, const unsigned int& numberOfPeriods);

int correctForDeadtimeBob(TH1D*& TOFtoCorrect, TH1D*& correctedTOF, const double deadtimeBins, const double deadtimeTransitionBins, const unsigned int& numberOfPeriods);

#endif /* CORRECT_FOR_DEADTIME_H */
