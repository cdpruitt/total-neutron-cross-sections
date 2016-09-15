#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

#include "TH1I.h"
#include <string>
#include <vector>

void scaleBins(std::vector<double> inputBins, int nInputBins, int scaledown, std::vector<double>& outputBins);

double tofToRKE(double TOF);

TH1* timeBinsToRKEBins(TH1I *inputHisto, std::string name);

#endif
