#ifndef FIT_DATA_H
#define FIT_DATA_H

class FitData
{
    double trigger1Time = 0;
    double trigger2Time = 0;
    double peak1Amplitude = 0;
    double peak2Amplitude = 0;
    double peak1Derivative = 0;
    double peak2Derivative = 0;
    double chiSquare = 0;
    bool goodFit = false;

    void clear()
    {
        trigger1Time = 0;
        trigger2Time = 0;
        peak1Amplitude = 0;
        peak2Amplitude = 0;
        peak1Derivative = 0;
        peak2Derivative = 0;
        chiSquare = 0;
        goodFit = false;
    }
}

#endif
