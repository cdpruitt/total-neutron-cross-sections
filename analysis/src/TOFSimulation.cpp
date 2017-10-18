#include <string>
#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TRandom3.h"

using namespace std;

const double DEADTIME_PERIOD = 10;
const double DEADTIME_TRANSITION_PERIOD = 50;

void incrementBin(unsigned int& j, unsigned int& deadtimeCounter)
{
    j++;

    if(deadtimeCounter>0)
    {
        deadtimeCounter--;
    }
}

int main(int argc, char** argv)
{
    string inputFileName = argv[1];
    string rateDistributionHistoName = argv[2];
    string outputFileName = argv[3];

    unsigned int numberOfPeriods = 1000000;

    TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
    if(!inputFile->IsOpen())
    {
        cerr << "Error: failed to open " << inputFileName << "  to fill histos." << endl;
        return 1;
    }
    
    TH1D* rateDistributionHisto = (TH1D*)inputFile->Get(rateDistributionHistoName.c_str());
    if(!rateDistributionHisto)
    {
        cerr << "Error: failed to open " << rateDistributionHistoName << " in " << inputFileName << " to read input rate distribution." << endl;
        return 1;
    }

    TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");

    unsigned int numberOfBins = rateDistributionHisto->GetNbinsX();

    vector<double> eventRate(numberOfBins);

    for(unsigned int i=0; i<numberOfBins; i++)
    {
        eventRate[i] = rateDistributionHisto->GetBinContent(i);
    }

    TRandom3* rng = new TRandom3();

    vector<unsigned int> currentPeriodEvents;
    vector<unsigned int> currentPeriodLiveEvents;

    TH1D* allEvents = new TH1D("allEvents","allEvents",numberOfBins,0,numberOfBins);
    TH1D* allLiveEvents = new TH1D("allLiveEvents","allEvents",numberOfBins,0,numberOfBins);

    unsigned int deadtimeCounter = 0;

    for(unsigned int i=0; i<numberOfPeriods; i++)
    {
        currentPeriodEvents.clear();
        currentPeriodLiveEvents.clear();

        // populate current period with events
        for(unsigned int j=0; j<numberOfBins; incrementBin(j, deadtimeCounter))
        {
            if(eventRate[j]>rng->Uniform(0,1))
            {
                // add an event to this bin for this period
                currentPeriodEvents.push_back(j);

                if(deadtimeCounter==0)
                {
                    // the detector is "live"
                    currentPeriodLiveEvents.push_back(j);
                }

                deadtimeCounter=DEADTIME_PERIOD;
            }
        }

        for(unsigned int b : currentPeriodEvents)
        {
            allEvents->Fill(b);
        }

        for(unsigned int b : currentPeriodLiveEvents)
        {
            allLiveEvents->Fill(b);
        }
    }

    allEvents->Write();
    allLiveEvents->Write();

    inputFile->Close();
    outputFile->Close();

    return 0;
}
