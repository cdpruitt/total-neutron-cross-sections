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

using namespace std;

extern RawEvent rawEvent;

// Header fields for events
unsigned long size;
unsigned long evtype;
unsigned long channel;
unsigned long timetag, timetagP, timeCoin;
unsigned int channelCoin = 1; // used to keep track of coincidence between channels

// Create histograms for DPP data
TH2S* FTLR;

double leftChannelTime;
double rightChannelTime;
double leftChannelFinetime;
double rightChannelFinetime;

void compareFineTimes(int numberOfEvents, TTree* tree)
{
    FTLR = new TH2S("outFTLR","outFTLR",1023,0,2,1023,0,2); // a ROOT plot to show fine time (FT) of events
    FTLR->SetMarkerStyle(20);

    // create text output to examine time differences from one event to the next
    ofstream timeDiff ("output/timeDiff.txt");

    int totalEntries = tree->GetEntries();
    if(totalEntries>numberOfEvents)
    {
        totalEntries = numberOfEvents;
    }

    for(int i=0; i<totalEntries; i++)
    {
        tree->GetEntry(i);
        if(rawEvent.chNo==6)
        {
            leftChannelTime = rawEvent.timetag;
            leftChannelFinetime = ((double)2000/1024)*rawEvent.fineTime;
        }

        if(rawEvent.chNo==7)
        {
            rightChannelTime = rawEvent.timetag;
            rightChannelFinetime = ((double)2000/1024)*rawEvent.fineTime;
        }

        if(leftChannelTime == rightChannelTime)
        {
            FTLR->Fill(leftChannelFinetime,rightChannelFinetime);
        }

        cout << "left finetime = " << rawEvent.chNo << endl;
    }
}

int main(int, char* argv[])
{
    cout << setprecision(10);

    int numberOfEvents = atoi(argv[2]);

    string rawFileName = argv[1];

    TFile* rawFile = new TFile(rawFileName.c_str(),"READ");
    TTree* tree = (TTree*)rawFile->Get("tree");

    if(!rawFile->Get("tree"))
    {
        cerr << "Error: failed to find raw tree in " << rawFileName << endl;
        exit(1);
    }

    // link the tree from the input file to our event variables
    tree->SetBranchAddress("chNo",&rawEvent.chNo);
    tree->SetBranchAddress("timetag",&rawEvent.timetag);
    tree->SetBranchAddress("fineTime",&rawEvent.fineTime);

    // create a ROOT file
    TFile *file; 
    file = new TFile("output/finetimeCheck.root","RECREATE");

    compareFineTimes(numberOfEvents, tree);

    file->Write();
}
