#include <iostream>
#include <string>
#include <vector>

#include "TH1.h"

#include "../include/config.h"
#include "../include/correctForDeadtime.h"

using namespace std;

// experimentally-determined digitizer deadtime
const int DEADTIME_PERIOD = 150; // in ns
const int DEADTIME_TRANSITION_PERIOD = 15; // in ns

double applyDeadtimeCorrection(TH1D* rawTOF, TH1D* correctedTOF, const vector<double>& deadtimesPerBin)
{
    if(deadtimesPerBin.size()==0)
    {
        cerr << "Error: cannot apply deadtime correction using empty deadtime correction list" << endl;
    }

    double sumOfDeadtimes = 0; // for computing average deadtime per bin

    for(int i=0; (size_t)i<deadtimesPerBin.size(); i++)
    {
        if(deadtimesPerBin[i]>=1)
        {
            cerr << "Error: attempted to correct for deadtime, but encountered deadtime >100% (deadtime was "
                << deadtimesPerBin[i] << "). Exiting." << endl;
            exit(1);
        }

        correctedTOF->SetBinContent(i+1,rawTOF->GetBinContent(i+1)/(1-deadtimesPerBin[i]));
        sumOfDeadtimes += deadtimesPerBin[i];
    }
    
    return sumOfDeadtimes/deadtimesPerBin.size();
}

int generateDeadtimeCorrection(TH1D* tof, unsigned int numberOfMacros,
        vector<double>& deadtimeCorrectionList)
{
    if(!tof)
    {
        cerr << "Error: empty TOF pointer in generateDeadtimeCorrection." << endl;
        return 1;
    }

    if(numberOfMacros==0)
    {
        cerr << "Error: number of macros was 0 while trying to "
            << "generate deadtime correction." << endl;
        return 1;
    }

    unsigned long int numberOfMicros = numberOfMacros
        *(config.facilityConfig.MACRO_LENGTH
        /(double)config.facilityConfig.MICRO_LENGTH);

    vector<double> eventsPerMicroPerBin(config.plotConfig.TOF_BINS);

    string name = tof->GetName();
    
    const int deadtimeBins = ((double)config.plotConfig.TOF_BINS/config.plotConfig.TOF_RANGE)*DEADTIME_PERIOD;

    const int deadtimeTransitionBins = ((double)config.plotConfig.TOF_BINS/config.plotConfig.TOF_RANGE)*DEADTIME_TRANSITION_PERIOD;

    for(int j=0; j<eventsPerMicroPerBin.size(); j++)
    {
        eventsPerMicroPerBin[j] = tof->GetBinContent(j+1)/((double)numberOfMicros);
    }

    // find the fraction of the time that the detector is dead for each bin in the micropulse
    // set deadtime fraction base case

    // use deadtime base case to calculate deadtime for remaining bins
    deadtimeCorrectionList.resize(eventsPerMicroPerBin.size());

    for(int j=0; j<config.plotConfig.TOF_BINS; j++)
    {
        for(int k=j-(deadtimeBins+deadtimeTransitionBins); k<j; k++)
        {
            if(k<(j-deadtimeBins))
            {
                if(k<0)
                {
                    deadtimeCorrectionList[j] += eventsPerMicroPerBin[k+config.plotConfig.TOF_BINS]*(double)(1-deadtimeCorrectionList[j])*((deadtimeBins+deadtimeTransitionBins-(j-k))/(double)deadtimeTransitionBins);
                }

                else
                {
                    deadtimeCorrectionList[j] += eventsPerMicroPerBin[k]*(double)(1-deadtimeCorrectionList[j])*((deadtimeBins+deadtimeTransitionBins-(j-k))/(double)deadtimeTransitionBins);
                }
            }

            else
            {
                if(k<0)
                {
                    deadtimeCorrectionList[j] += eventsPerMicroPerBin[k+config.plotConfig.TOF_BINS]*(double)(1-deadtimeCorrectionList[j]);
                }

                else
                {
                    deadtimeCorrectionList[j] += eventsPerMicroPerBin[k]*(double)(1-deadtimeCorrectionList[j]);
                }
            }
        }

        deadtimeCorrectionList[j] += (eventsPerMicroPerBin[j]/2)*(1-deadtimeCorrectionList[j]); // last bin contributes 1/2 its value
    }

    return 0;
}
/*

// perform iterative deadtime correction, until average deadtime changes
    // <0.1%
    vector<double> deadtimeBins = generateDeadtimeCorrection(TOFHistos[i], microsPerTarget[i]);

        double averageDeadtimeDiff = applyDeadtimeCorrection(TOFHistos[i], deadtimeBins) - 0;

        while(averageDeadtimeDiff>0.001)
          {
          rawTOF = correctedTOF;

          vector<double> prevDeadtimeBins = deadtimeBins;
          deadtimeBins = generateDeadtimeCorrection(rawTOF, microsPerTarget[i]);

          for(unsigned int j=0; j<deadtimeBins.size(); j++)
          {
          deadtimeBins[j] = deadtimeBins[j]-prevDeadtimeBins[j];
          }

          correctedTOF = ((TH1D*)plots[i]->getTOFHisto());
          deadtimeH = ((TH1D*)plots[i]->getDeadtimeHisto());
          averageDeadtimeDiff = applyDeadtimeCorrection(rawTOF, correctedTOF, deadtimeH, deadtimeBins) - averageDeadtimeDiff;
          }
    }

    // calculate likelihood of double peak in wavelet, for each TOF bin
    //const int waveletBins = (config.plotConfig.TOF_BINS/config.plotConfig.TOF_RANGE)*WAVELET_PERIOD;

*/

