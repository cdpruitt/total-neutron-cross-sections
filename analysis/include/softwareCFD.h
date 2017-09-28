#ifndef SOFTWARE_CFD_H
#define SOFTWARE_CFD_H

double calculateCFDTime(std::vector<int>* waveform,
        double baseline,
        double fraction,
        unsigned int delay,
        bool isPositiveSignal);

#endif /* SOFTWARE_CFD_H */
