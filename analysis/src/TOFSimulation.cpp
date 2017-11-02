#include <string>
#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TRandom2.h"

#include "../include/correctForDeadtime.h"
#include "../include/plots.h"

using namespace std;

const double DEADTIME_BINS = 100;
const double DEADTIME_TRANSITION_BINS = 0;
const unsigned int PERIOD_RESET_NUMBER = 250;

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
    unsigned int numberOfPeriods = atoi(argv[3]);
    string outputFileName = argv[4];

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
        eventRate[i] = rateDistributionHisto->GetBinContent(i+1);
    }

    TRandom2* rng = new TRandom2();

    vector<unsigned int> currentPeriodEvents;
    vector<unsigned int> currentPeriodLiveEvents;

    TH1D* allEvents = new TH1D("allEvents","allEvents",numberOfBins,0,numberOfBins);
    TH1D* allLiveEvents = new TH1D("allLiveEvents","allEvents",numberOfBins,0,numberOfBins);

    unsigned int deadtimeCounter = 0;

    unsigned int eventsToAdd = 0;

    for(unsigned int i=0; i<numberOfPeriods; i++)
    {
        currentPeriodEvents.clear();
        currentPeriodLiveEvents.clear();

        // populate current period with events
        for(unsigned int j=0; j<numberOfBins; incrementBin(j, deadtimeCounter))
        {
            eventsToAdd = rng->Poisson(eventRate[j]);

            if(eventsToAdd>0)
            {
                // add events to this bin for this period
                while(eventsToAdd>0)
                {
                    currentPeriodEvents.push_back(j);
                    eventsToAdd--;
                }

                if(deadtimeCounter==0)
                {
                    // the detector is "live"
                    currentPeriodLiveEvents.push_back(j);

                    deadtimeCounter=DEADTIME_BINS;
                }
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

        if(i%1000==0)
        {
            cout << "Ran simulation through " << i << " number of periods.\r";
            fflush(stdout);
        }

        if(i%PERIOD_RESET_NUMBER==0)
        {
            // reached end of "macropulse"; reset the deadtime counter for the
            // next "macropulse" start
            deadtimeCounter = 0;
        }
    }

    TH1D* measured = (TH1D*)allLiveEvents->Clone("measured");
    TH1D* corrected = (TH1D*)allLiveEvents->Clone("corrected");

    correctForDeadtimeBob(measured, corrected, DEADTIME_BINS, DEADTIME_TRANSITION_BINS, numberOfPeriods);

    corrected->Write();

    allEvents->Write();
    allLiveEvents->Write();

    inputFile->Close();
    outputFile->Close();

    return 0;
}
