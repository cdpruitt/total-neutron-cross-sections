#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include "TH1S.h"
#include "TH2S.h"
#include "TFile.h"
#include "TTree.h"

#include "../include/dataStructures.h"
#include "../include/raw.h"
#include "../include/experiment.h"

using namespace std;

const unsigned int TIME_CHECK_TOLERANCE = 10; // in ns
const unsigned int NUMBER_OF_EVENTS = 10000;  // number of events to examine for time correlation

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

    // create a ROOT file for holding time correlation plots
    TFile* file = new TFile("timeCheckOutput/detTimeCheck.root","RECREATE");

    TH2D* detTimeCorrelation = new TH2D("left/right time correlation","left/right time correlation",1000,43,53,1000,43,53);
    detTimeCorrelation->GetXaxis()->SetTitle("right detector");
    detTimeCorrelation->GetYaxis()->SetTitle("left detector");
    detTimeCorrelation->SetMarkerStyle(20);

    TH1D* detTimeDifference = new TH1D("time difference","diffLR",1000,-3,3);
    detTimeDifference->GetXaxis()->SetTitle("left detector time - right detector time");
    detTimeDifference->SetMarkerStyle(20);

    unsigned long numberOfEventsProcessed = 0;

    double ch6Timetag = 0;
    double ch6FineTime = 0;

    double ch7Timetag = 0;
    double ch7FineTime = 0;

    RawEvent rawEvent;

    while(!inFile.eof() && numberOfEventsProcessed < NUMBER_OF_EVENTS)
    {
        if(readEvent(inFile, rawEvent))
        {
            if(rawEvent.chNo==6)
            {
                ch6Timetag = rawEvent.timetag;
                ch6FineTime = calculateCFDTime(rawEvent.waveform,
                            rawEvent.baseline,
                            CFD_FRACTION,
                            CFD_DELAY,
                            false);
            }

            else if(rawEvent.chNo==7)
            {
                ch7Timetag = rawEvent.timetag;
                ch7FineTime = calculateCFDTime(rawEvent.waveform,
                            rawEvent.baseline,
                            CFD_FRACTION,
                            CFD_DELAY,
                            false);
            }
        }

        if(abs(ch6Timetag-ch7Timetag)<=TIME_CHECK_TOLERANCE && ch6FineTime>=0 && ch7FineTime>=0)
        {
            detTimeCorrelation->Fill(ch7FineTime, ch6FineTime+(ch6Timetag-ch7Timetag));
            detTimeDifference->Fill((ch6FineTime+ch6Timetag)-(ch7FineTime+ch7Timetag));
        }
    }

    file->Write();

    return 0;
}
