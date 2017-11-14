#ifndef CORRECT_FOR_DEADTIME_H
#define CORRECT_FOR_DEADTIME_H

#include <vector>

#include "TH1.h"
#include "TDirectory.h"

int correctForDeadtime(TH1D*& TOFtoCorrect, TH1D*& correctedTOF, const unsigned int& numberOfPeriods);

#endif /* CORRECT_FOR_DEADTIME_H */
