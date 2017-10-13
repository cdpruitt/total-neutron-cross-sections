#ifndef SOFTWARE_CFD_H
#define SOFTWARE_CFD_H

#include <vector>

double calculateCFDTime(const std::vector<int>& waveform,
        const double& baseline,
        const double& fraction,
        const unsigned int& delay);

double calculateMacropulseFineTime(std::vector<int>* waveform, double threshold);

#endif /* SOFTWARE_CFD_H */
