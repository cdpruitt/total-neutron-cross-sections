#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

#include "TH1S.h"
#include "TH2S.h"
#include "TFile.h"
#include "TTree.h"

#include "../include/dataStructures.h"
#include "../include/raw.h"
#include "../include/experiment.h"
#include "../include/experimentalConfig.h"
#include "../include/softwareCFD.h"

using namespace std;

const unsigned int TIME_CHECK_TOLERANCE = 5; // in ns
const unsigned int NUMBER_OF_EVENTS = 1000000;  // number of events to examine for time correlation

ExperimentalConfig experimentalConfig;

int main(int argc, char** argv)
{
    // open input file
    string inFileName = argv[1];
    ifstream inFile;
    inFile.open(inFileName,ios::binary);

    if (!inFile)
    {
        cout << "Failed to open " << inFileName << ". Please check that the file exists" << endl;
        exit(1);
    }

    cout << inFileName << " opened successfully." << endl;

    string outFileLocation = argv[2];
    outFileLocation = outFileLocation + "detTimeCheck.root";

    // create a ROOT file for holding time correlation plots
    TFile* outFile = new TFile(outFileLocation.c_str(),"RECREATE");

    TH2D* detTimeCorrelation = new TH2D("LR time correlation","LR time correlation",300,30,60,300,30,60);
    detTimeCorrelation->GetXaxis()->SetTitle("right detector");
    detTimeCorrelation->GetYaxis()->SetTitle("left detector");
    detTimeCorrelation->SetMarkerStyle(7);

    TH1D* detTimeDifference = new TH1D("LR time difference","LR time difference",1000,-5,5);
    detTimeDifference->GetXaxis()->SetTitle("left detector time - right detector time");
    detTimeDifference->SetMarkerStyle(7);

    TH1D* ch6FineTimeH = new TH1D("ch6 fine time", "ch6 fine time", 600, 0, 60);
    TH1D* ch7FineTimeH = new TH1D("ch7 fine time", "ch7 fine time", 600, 0, 60);

    unsigned long numberOfEventsProcessed = 0;

    double ch6Timetag = 0;
    double ch6FineTime = 0;

    double ch7Timetag = 0;
    double ch7FineTime = 0;

    unsigned long numberOfBadCh6FineTime = 0;
    unsigned long numberOfBadCh7FineTime = 0;
    unsigned long numberOfCh6Events = 0;
    unsigned long numberOfCh7Events = 0;

    RawEvent rawEvent;

    while(!inFile.eof() && numberOfEventsProcessed < NUMBER_OF_EVENTS)
    {
        if(readEvent(inFile, rawEvent))
        {
            if(rawEvent.chNo==6)
            {
                ch6Timetag = rawEvent.timetag*experimentalConfig.timeConfig.SAMPLE_PERIOD;
                ch6FineTime = calculateCFDTime(
                        rawEvent.waveform,
                        rawEvent.baseline,
                        experimentalConfig.timeConfig.CFD_FRACTION,
                        experimentalConfig.timeConfig.CFD_DELAY,
                        false
                        )*experimentalConfig.timeConfig.SAMPLE_PERIOD;

                ch6FineTimeH->Fill(ch6FineTime);

                if(ch6FineTime<0)
                {
                    stringstream name;
                    name << "Event" << numberOfEventsProcessed;
                    //TH1D* badFineTimeHisto = new TH1D(name.str().c_str(),name.str().c_str(),rawEvent.waveform->size(),0,rawEvent.waveform->size());
                    /*for(unsigned int i=0; i<rawEvent.waveform->size(); i++)
                    {
                        badFineTimeHisto->SetBinContent(i+1,rawEvent.waveform->at(i));
                    }
                    badFineTimeHisto->Write();
                    */

                    numberOfBadCh6FineTime++;
                }

                numberOfCh6Events++;
            }

            else if(rawEvent.chNo==7)
            {
                ch7Timetag = rawEvent.timetag*experimentalConfig.timeConfig.SAMPLE_PERIOD;
                ch7FineTime = calculateCFDTime(
                        rawEvent.waveform,
                        rawEvent.baseline,
                        experimentalConfig.timeConfig.CFD_FRACTION,
                        experimentalConfig.timeConfig.CFD_DELAY,
                        false
                        )*experimentalConfig.timeConfig.SAMPLE_PERIOD;

                ch7FineTimeH->Fill(ch7FineTime);

                if(ch7FineTime<0)
                {
                    numberOfBadCh7FineTime++;
                }

                numberOfCh7Events++;
            }
        }

        if(abs(ch6Timetag-ch7Timetag)<=TIME_CHECK_TOLERANCE && ch6FineTime>=0 && ch7FineTime>=0)
        {
            detTimeCorrelation->Fill(ch7FineTime, ch6FineTime+(ch6Timetag-ch7Timetag));
            detTimeDifference->Fill((ch6FineTime+ch6Timetag)-(ch7FineTime+ch7Timetag));
        }

        if(numberOfEventsProcessed%10000==0)
        {
            cout << "Processed " << numberOfEventsProcessed << " events through time check\r";
            fflush(stdout);
        }

        numberOfEventsProcessed++;
    }

    cout << endl;
    cout << "Finished processing " << numberOfEventsProcessed << " events through time check." << endl;

    cout << "Successfully calculated a fine time on " << 100*(double)(numberOfCh6Events-numberOfBadCh6FineTime)/(numberOfCh6Events) << "% of ch6 events." << endl;
    cout << "Successfully calculated a fine time on " << 100*(double)(numberOfCh7Events-numberOfBadCh7FineTime)/(numberOfCh7Events) << "% of ch7 events." << endl;

    outFile->Write();
    outFile->Close();

    return 0;
}
