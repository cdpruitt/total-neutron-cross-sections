#ifndef SOFTWARE_CFD_H
#define SOFTWARE_CFD_H

#include <vector>

double calculateCFDTime(std::vector<int>* waveform,
        double baseline,
        double fraction,
        unsigned int delay);

double calculateMacropulseFineTime(std::vector<int>* waveform, double threshold);

#endif /* SOFTWARE_CFD_H */
