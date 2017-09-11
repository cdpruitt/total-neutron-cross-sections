// Reads output from washudaq-X.X and sorts into ROOT and text files
// The expected event file structure is as follows:
//
//      EVENT HEADER:
//
//      Size            |   uint32; number of bytes in the event (self-inclusive)
//      Event Type      |   uint32; possible values are 1 (DPP data) or 2 (waveform)
//      Channel         |   uint32; channel number of event, starting with 0
//      Time Tag        |   uint32; coarse time (2 ns units) of event trigger
//
//      (every event has an EVENT HEADER)
//
//      ---------------------------------------------------------------------------
//
//      DPP EVENT BODY:
//
//      Zero Crossing   |   uint16; interpolated zero-crossing, in picoseconds
//      Short Gate Q    |   uint16; charge integrated in the short gate
//      Long Gate Q     |   uint16; charge integrated in the long gate
//      Baseline        |   uint16; baseline level of ADC
//      Pile up rej.    |   uint16; pile-up rejection flag (not yet implemented)
//      Num. of samples |   uint32; number of waveform samples to follow (mixed mode)
//      Samples         |   series of uint16 giving raw waveform samples
//
//      (events of type 1 have DPP EVENT BODY)
//
//      ---------------------------------------------------------------------------
//
//      WAVEFORM EVENT BODY
//
//      Num. of samples |   uint32; number of waveform samples to follow (mixed mode)
//      Samples         |   series of uint16 giving raw waveform samples
//
//      (events of type 2 have WAVEFORM EVENT BODY)

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include "TH1S.h"
#include "TH2S.h"
#include "TFile.h"
#include "TTree.h"

#include "../include/dataStructures.h"
#include "../include/raw.h"
#include "../include/physicalConstants.h"
#include "../include/analysisConstants.h"
#include "../include/experiment.h"
#include "../include/digitizerConstants.h"


using namespace std;

extern RawEvent rawEvent;

// Header fields for events
unsigned long size;
unsigned long evtype;
unsigned long channel;
unsigned long timetag, timetagP, timeCoin;
unsigned int channelCoin = 1; // used to keep track of coincidence between channels

// Create histograms for DPP data
TH2I* FTLR;
TH1I* diffLR;

const unsigned int NUMBER_OF_EVENTS = 10000;

int main(int argc, char** argv)
{
    string inFileName = argv[1];
    ifstream inFile;
    inFile.open(inFileName,ios::binary);

    if (!inFile)
    {
        cout << "Failed to open " << inFileName << ". Please check that the file exists" << endl;
        exit(1);
    }

    cout << inFileName << " opened successfully." << endl;

    // create a ROOT file for holding time check plots
    TFile *file; 
    file = new TFile("timeCheckOutput/fineTimeCheck.root","RECREATE");

    FTLR = new TH2I("leftDet","rightDet",1000,43,53,1000,43,53); // a ROOT plot to show fine time (FT) of events
    FTLR->SetMarkerStyle(20);

    diffLR = new TH1I("diffLR","diffLR",1000,-3,3); // a ROOT plot to show fine time (FT) of events
    FTLR->SetMarkerStyle(20);


    // create text output to examine time differences from one event to the next
    ofstream timeDiff ("timeCheckOutput/timeDiff.txt");

    unsigned int numberOfEventsProcessed = 0;

    unsigned int ch6Timetag = 0;
    double ch6FineTime = 0;

    unsigned int ch7Timetag = 0;
    double ch7FineTime = 0;

    while(!inFile.eof() && numberOfEventsProcessed < NUMBER_OF_EVENTS)
    {
        if(readEvent(inFile))
        {
            if(rawEvent.chNo==6)
            {
                ch6Timetag = rawEvent.timetag;
                ch6FineTime = calculateCFDTime(rawEvent.waveform,
                            rawEvent.baseline,
                            CFD_FRACTION,
                            CFD_DELAY,
                            false)-DETECTOR_PRETRIGGER_SAMPLES;
            }

            if(rawEvent.chNo==7)
            {
                ch7Timetag = rawEvent.timetag;
                ch7FineTime = calculateCFDTime(rawEvent.waveform,
                            rawEvent.baseline,
                            CFD_FRACTION,
                            CFD_DELAY,
                            false)-DETECTOR_PRETRIGGER_SAMPLES;
            }

            if(abs((int)ch6Timetag-(int)ch7Timetag)<=TIME_CHECK_TOLERANCE && ch6FineTime>=0 && ch7FineTime>=0)
            {
                FTLR->Fill(ch6FineTime, ch7FineTime+(ch7Timetag-ch6Timetag));
                diffLR->Fill((ch7FineTime+ch7Timetag)-(ch6FineTime+ch6Timetag));
            }
        }
    }

    file->Write();

    return 0;
}
