#include <iostream>
#include <string>
#include <vector>

#include "TH1.h"
#include "TDirectory.h"

#include "../include/correctForDeadtime.h"
#include "../include/config.h"

using namespace std;

extern Config config;

// logistic curve for deadtime response
// (calculated by fitting "time difference between events" histogram using a
// logistic curve)
const double DEADTIME_LOGISTIC_k = 0.546698;
const double DEADTIME_LOGISTIC_MU = 159.73; // in ns

// input x in ns
double logisticDeadtimeFunction(double x)
{
    if(x<140)
    {
        return 1;
    }

    return 1-1/(1+exp(-DEADTIME_LOGISTIC_k*(x-DEADTIME_LOGISTIC_MU)));
}

int correctForDeadtime(TH1D*& TOFtoCorrect, TH1D*& correctedTOF, const unsigned int& numberOfPeriods)
{
    int DEADTIME_BINS = 170*config.plot.TOF_BINS_PER_NS;

    string targetName = TOFtoCorrect->GetName();
    unsigned int numberOfBins = TOFtoCorrect->GetNbinsX();

    vector<double> measuredRatePerBin(numberOfBins);

    for(int i=0; i<measuredRatePerBin.size(); i++)
    {
        measuredRatePerBin[i] = TOFtoCorrect->GetBinContent(i+1)/(double)numberOfPeriods;

        /*if(i>85*config.plot.TOF_BINS_PER_NS && i<95*config.plot.TOF_BINS_PER_NS)
        {
            measuredRatePerBin[i] *= 1.2;
        }*/
    }

    vector<long double> trueRatePerBin(numberOfBins);

    for(int i=0; i<numberOfBins; i++)
    {
        double cumulativeDeadtime = 0;
        for(int j=0; j<DEADTIME_BINS; j++)
        {
            if((i-j)<0)
            {
                cumulativeDeadtime += measuredRatePerBin[(i-j)+numberOfBins]*logisticDeadtimeFunction(((double)j)/config.plot.TOF_BINS_PER_NS);
            }

            else
            {
                cumulativeDeadtime += measuredRatePerBin[i-j]*logisticDeadtimeFunction(((double)j)/config.plot.TOF_BINS_PER_NS);
            }
        }

        trueRatePerBin[i] = -log(1-(measuredRatePerBin[i])/(1-cumulativeDeadtime));
    }

    for(int i=0; i<numberOfBins; i++)
    {
        correctedTOF->SetBinContent(i+1, trueRatePerBin[i]*numberOfPeriods);
    }

    cout << "Finished deadtime correction for " << targetName << "." << endl;

    return 0;
}
