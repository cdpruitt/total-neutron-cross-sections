#include <iostream>
#include <string>
#include <vector>

#include "TH1.h"
#include "TDirectory.h"

#include "../include/correctForDeadtime.h"

using namespace std;

const double CONVERGENCE_CUTOFF = 0.99999;

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

int generateDeadtimeCorrection(TH1D* tof, double deadtimeBins, double deadtimeTransitionBins, unsigned int numberOfPeriods, vector<double>& deadtimeCorrection)
{
    cout << "Generating deadtime correction for target " << tof->GetName() << "\r";

    if(!tof)
    {
        cerr << "Error: empty TOF pointer in generateDeadtimeCorrection." << endl;
        return 1;
    }

    if(numberOfPeriods==0)
    {
        cerr << "Error: number of macros was 0 while trying to "
            << "generate deadtime correction." << endl;
        return 1;
    }

    string name = tof->GetName();
    unsigned int numberOfBins = tof->GetNbinsX();

    vector<double> eventsPerMicroPerBin(numberOfBins);

    for(int j=0; j<eventsPerMicroPerBin.size(); j++)
    {
        eventsPerMicroPerBin[j] = tof->GetBinContent(j+1)/numberOfPeriods;
    }

    // find the fraction of the time that the detector is dead for each bin in the micropulse
    // set deadtime fraction base case

    // use deadtime base case to calculate deadtime for remaining bins
    deadtimeCorrection.resize(eventsPerMicroPerBin.size());

    for(int j=0; j<numberOfBins; j++)
    {
        for(int k=j-(deadtimeBins+deadtimeTransitionBins); k<j; k++)
        {
            if(k<(j-deadtimeBins))
            {
                if(k<0)
                {
                    deadtimeCorrection[j] += eventsPerMicroPerBin[k+numberOfBins]*(double)(1-deadtimeCorrection[j])*((deadtimeBins+deadtimeTransitionBins-(j-k))/(double)deadtimeTransitionBins);
                }

                else
                {
                    deadtimeCorrection[j] += eventsPerMicroPerBin[k]*(double)(1-deadtimeCorrection[j])*((deadtimeBins+deadtimeTransitionBins-(j-k))/(double)deadtimeTransitionBins);
                }
            }

            else
            {
                if(k<0)
                {
                    deadtimeCorrection[j] += eventsPerMicroPerBin[k+numberOfBins]*(double)(1-deadtimeCorrection[j]);
                }

                else
                {
                    deadtimeCorrection[j] += eventsPerMicroPerBin[k]*(double)(1-deadtimeCorrection[j]);
                }
            }
        }

        deadtimeCorrection[j] += (eventsPerMicroPerBin[j]/2)*(1-deadtimeCorrection[j]); // last bin contributes 1/2 its value
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

int correctForDeadtime(TH1D*& rawTOF, TH1D*& correctedTOF, const double deadtimeBins, const double deadtimeTransitionBins, const unsigned int& numberOfMacros)
{
    string targetName = rawTOF->GetName();

    unsigned int iteration = 0;

    rawTOF->Reset();

    while(!testDeadtimeConvergence(rawTOF, correctedTOF))
    {
        iteration++;

        TH1D* inputTOF = (TH1D*)correctedTOF->Clone();
        inputTOF->Add(rawTOF,-1);

        vector<double> deadtimeCorrection;

        generateDeadtimeCorrection(inputTOF, deadtimeBins, deadtimeTransitionBins, numberOfMacros, deadtimeCorrection);

        string deadtimeName = targetName + "Deadtime" + to_string(iteration);

        unsigned int numberOfBins = inputTOF->GetNbinsX();
        TH1D* deadtimeHisto = new TH1D(deadtimeName.c_str(), deadtimeName.c_str(),numberOfBins,0,numberOfBins);

        for(int j=0; j<deadtimeCorrection.size(); j++)
        {
            deadtimeHisto->SetBinContent(j+1,deadtimeCorrection[j]);
        }

        deadtimeHisto->Write();

        rawTOF = correctedTOF;

        string correctedName = targetName + "Iteration" + to_string(iteration);
        correctedTOF = (TH1D*)correctedTOF->Clone(correctedName.c_str());

        applyDeadtimeCorrection(correctedTOF, deadtimeCorrection);

        correctedTOF->Write();
    }

    cout << endl;
    cout << "Finished deadtime correction for " << targetName << "TOF." << endl;

    return 0;
}

int correctForDeadtime2(TH1D*& TOFtoCorrect, TH1D*& correctedTOF, const double deadtimeBins, const double deadtimeTransitionBins, const unsigned int& numberOfPeriods)
{
    string targetName = TOFtoCorrect->GetName();
    unsigned int numberOfBins = TOFtoCorrect->GetNbinsX();

    vector<double> measuredRatePerBin(numberOfBins);

    for(int i=0; i<measuredRatePerBin.size(); i++)
    {
        measuredRatePerBin[i] = TOFtoCorrect->GetBinContent(i+1)/(double)numberOfPeriods;
    }

    vector<long double> trueRatePerBin(numberOfBins);
    vector<long double> deadtimeCorrection(numberOfBins);

    for(int i=0; i<numberOfBins; i++)
    {
        for(int j=i-(deadtimeBins+deadtimeTransitionBins); j<i; j++)
        {
            int k=i-(deadtimeBins+deadtimeTransitionBins);

            if((j-k)<deadtimeTransitionBins)
            {
                // deadtime transition region
                if(j<0)
                {
                    deadtimeCorrection[i] += trueRatePerBin[j+numberOfBins]*(double)(1-deadtimeCorrection[i])*((j-k)/(double)deadtimeTransitionBins);
                }

                else
                {
                    deadtimeCorrection[i] += trueRatePerBin[j]*(double)(1-deadtimeCorrection[i])*((j-k)/(double)deadtimeTransitionBins);
                }
            }

            else
            {
                if(j<0)
                {
                    deadtimeCorrection[i] += trueRatePerBin[j+numberOfBins]*(double)(1-deadtimeCorrection[i]);
                }

                else
                {
                    deadtimeCorrection[i] += trueRatePerBin[j]*(double)(1-deadtimeCorrection[i]);
                }
            }

            // calculate true event rate in this bin (current bin contributes half
            // its value to the deadtime of itself)
            if(deadtimeCorrection[i]>=1)
            {
                cerr << "Error: deadtime correction was calculated to be >1. Ending deadtime correction calculation for bin " << i << "." << endl;
                break;
            }

            trueRatePerBin[i] = 1 - sqrt(1+(2*measuredRatePerBin[i])/(deadtimeCorrection[i]-1));
        }

        cout << "bin i = " << i << ", true rate = " << trueRatePerBin[i] << endl;
    }

    // write out calculated deadtime and true event rate
    string deadtimeName = targetName + "Deadtime";

    TH1D* deadtimeHisto = new TH1D(deadtimeName.c_str(), deadtimeName.c_str(),numberOfBins,0,numberOfBins);

    for(int j=0; j<deadtimeCorrection.size(); j++)
    {
        deadtimeHisto->SetBinContent(j+1,deadtimeCorrection[j]);
    }

    deadtimeHisto->Write();

    for(int i=0; i<numberOfBins; i++)
    {
        correctedTOF->SetBinContent(i+1, trueRatePerBin[i]*numberOfPeriods);
    }

    cout << endl;
    cout << "Finished deadtime correction for " << targetName << "TOF." << endl;

    return 0;
}

int correctForDeadtimeBob(TH1D*& TOFtoCorrect, TH1D*& correctedTOF, const double deadtimeBins, const double deadtimeTransitionBins, const unsigned int& numberOfPeriods)
{
    string targetName = TOFtoCorrect->GetName();
    unsigned int numberOfBins = TOFtoCorrect->GetNbinsX();

    vector<double> measuredRatePerBin(numberOfBins);

    for(int i=0; i<measuredRatePerBin.size(); i++)
    {
        measuredRatePerBin[i] = TOFtoCorrect->GetBinContent(i+1)/(double)numberOfPeriods;
    }

    vector<long double> trueRatePerBin(numberOfBins);

    for(int i=0; i<numberOfBins; i++)
    {
        double cumulativeDeadtime = 0;
        for(int j=0; j<deadtimeBins+deadtimeTransitionBins; j++)
        {
            if(j>deadtimeBins)
            {
                if((i-j)<0)
                {
                    cumulativeDeadtime += measuredRatePerBin[(i-j)+numberOfBins]
                        *((deadtimeTransitionBins-(j-deadtimeBins))
                                /(double)deadtimeTransitionBins);
                }

                else
                {
                    cumulativeDeadtime += measuredRatePerBin[i-j]
                        *((deadtimeTransitionBins-(j-deadtimeBins))
                                /(double)deadtimeTransitionBins);
                }
            }

            else
            {
                if((i-j)<0)
                {
                    cumulativeDeadtime += measuredRatePerBin[(i-j)+numberOfBins];
                }

                else
                {
                    cumulativeDeadtime += measuredRatePerBin[i-j];
                }
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
