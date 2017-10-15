#include <iostream>
#include <string>
#include <vector>

#include "TH1.h"
#include "TDirectory.h"

#include "../include/config.h"
#include "../include/correctForDeadtime.h"

using namespace std;

// experimentally-determined digitizer deadtime
const int DEADTIME_PERIOD = 150; // in ns
const int DEADTIME_TRANSITION_PERIOD = 15; // in ns

const double CONVERGENCE_CUTOFF = 0.999;

bool testDeadtimeConvergence(TH1D* inputTOF, TH1D* correctedTOF)
{
    unsigned int numberOfBins = inputTOF->GetNbinsX();

    for(unsigned int i=1; i<=numberOfBins; i++)
    {
        // test that histos are within a convergence cutoff for every bin
        if((correctedTOF->GetBinContent(i) > 0) &&
           (inputTOF->GetBinContent(i)/correctedTOF->GetBinContent(i))<
           CONVERGENCE_CUTOFF)
        {
            return false;
        }
    }

    return true;
}

int generateDeadtimeCorrection(TH1D* tof, unsigned int numberOfMacros,
        vector<double>& deadtimeCorrectionList)
{
    cout << "Generating deadtime correction for target " << tof->GetName()<< endl;

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

void applyDeadtimeCorrection(TH1D* correctedTOF, const vector<double>& deadtimesPerBin)
{
    if(deadtimesPerBin.size()==0)
    {
        cerr << "Error: cannot apply deadtime correction using empty deadtime correction list" << endl;
        return;
    }

    for(int i=0; (size_t)i<deadtimesPerBin.size(); i++)
    {
        if(deadtimesPerBin[i]>=1)
        {
            cerr << "Error: attempted to correct for deadtime, but encountered deadtime >100% (deadtime was "
                << deadtimesPerBin[i] << "). Exiting." << endl;
            return;
        }

        correctedTOF->SetBinContent(i+1,correctedTOF->GetBinContent(i+1)/(1-deadtimesPerBin[i]));
    }

    return;
}

int correctForDeadtime(TH1D*& rawTOF, TH1D*& correctedTOF, const unsigned int& numberOfMacros, TDirectory*& outputDirectory)
{
    outputDirectory->cd();

    string targetName = rawTOF->GetName();

    unsigned int iteration = 0;

    rawTOF->Reset();

    while(!testDeadtimeConvergence(rawTOF, correctedTOF))
    {
        iteration++;

        TH1D* inputTOF = (TH1D*)correctedTOF->Clone();
        inputTOF->Add(rawTOF,-1);

        vector<double> deadtimeCorrectionList;

        generateDeadtimeCorrection(inputTOF, numberOfMacros, deadtimeCorrectionList);

        string deadtimeName = targetName + "Deadtime" + to_string(iteration);
        TH1D* deadtimeHisto = new TH1D(deadtimeName.c_str(), deadtimeName.c_str(),config.plotConfig.TOF_BINS,config.plotConfig.TOF_LOWER_BOUND,config.plotConfig.TOF_UPPER_BOUND);

        for(int j=0; j<deadtimeCorrectionList.size(); j++)
        {
            deadtimeHisto->SetBinContent(j+1,deadtimeCorrectionList[j]);
        }

        deadtimeHisto->Write();

        rawTOF = correctedTOF;

        string correctedName = targetName + "Iteration" + to_string(iteration);
        correctedTOF = (TH1D*)correctedTOF->Clone(correctedName.c_str());

        applyDeadtimeCorrection(correctedTOF, deadtimeCorrectionList);

        correctedTOF->Write();
    }

    return 0;
}

/*void applyDeadtimeCorrection(TH2D* rawTOF, TH2D* correctedTOF, const vector<double>& deadtimesPerBin)
  {
  unsigned int numberXBins = correctedTOF->GetNbinsX();
  unsigned int numberYBins = correctedTOF->GetNbinsY();

  if(numberXBins!=deadtimesPerBin.size() || numberYBins!=deadtimesPerBin.size())
  {
  cerr << "Error: cannot apply deadtime correction to histogram with different number of bins." << endl;
  return;
  }

  for(int i=0; i<numberXBins; i++)
  {
  for(int j=0; j<numberYBins; j++)
  {
  correctedTOF->SetBinContent(i+1,j+1,rawTOF->GetBinContent(i+1)/(1-deadtimesPerBin[i]));
  }
  }

  return;
  }*/


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

