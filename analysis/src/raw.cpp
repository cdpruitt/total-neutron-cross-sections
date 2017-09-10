/******************************************************************************
  raw.cpp
 ******************************************************************************/

// raw.cpp reads raw event data from a binary .evt file and outputs a ROOT file
// with a single ROOT tree containing all the raw event data.
// 
// The input .evt file is expected to consist of a series of events back-to-back,
// each with the following structure:
//
//   First, an EVENT HEADER:
//   ---------------------------------------------------------------------------
//   Size            |   uint32; number of bytes in the event (self-inclusive)
//   Event Type      |   uint32; possible values are 1 (DPP) or 2 (waveform)
//   Channel         |   uint32; channel number of event (range: 0-7)
//   Time Tag        |   uint32; coarse time (2 ns units) of event trigger
//   ---------------------------------------------------------------------------
//
//   Second, an EVENT BODY, either DPP (event type 1) or WAVEFORM (event type 2):
//   ---------------------------------------------------------------------------
//   DPP event body:
//
//   Extra select    |   uint16; describes the contents of Extras
//   Extras          |   uint32; additional PSD data (see description below)
//   Short Gate Q    |   uint16; charge integrated in the short gate
//   Long Gate Q     |   uint16; charge integrated in the long gate
//   Baseline        |   uint16; baseline level of ADC
//   Pile up rej.    |   uint16; pile-up rejection flag (not yet implemented)
//   Probe info      |   uint16; turns on an analog probe for examining waveform
//   Num. of samples |   uint32; number of waveform samples collected in Samples
//   Samples         |   series of uint16 giving raw waveform samples
//   ---------------------------------------------------------------------------
//   WAVEFORM event body:
//
//   Num. of samples |   uint32; number of waveform samples collected in Samples
//   Samples         |   series of uint16 giving raw waveform samples
//   ---------------------------------------------------------------------------
//
//
//   In the DPP event body, "Extras" is a 32-bit word whose content varies with
//   the value of EXTRA_SELECT as follows:
//   ---------------------------------------------------------------------------
//   [DEFAULT]  0: extended timestamp (bits 16-31) and baseline*4 (bits 0-15)
//              1: extended timestamp (bits 16-31) and flags (bits 0-15?)
//              2: extended timestamp (bits 16-31), flags (bits 10-15), and fine
//                 timestamp (bits 0-9)
//              3: pulse peak value (bits 0-15) [UNTESTED FEATURE]
//              5: CFD positive ZC (bits 16-31) and negative ZC (bits 0-15)
//              7: fixed value of 0x12345678
//   ---------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "TFile.h"
#include "TTree.h"

#include "../include/buffer.h"
#include "../include/raw.h"
#include "../include/dataStructures.h"
#include "../include/physicalConstants.h"
#include "../include/digitizerConstants.h"
#include "../include/branches.h"

using namespace std;

extern RawEvent rawEvent; // struct for storing raw event data;
                          // defined in branches.cpp

// read a word off the input file and store in a variable
bool readWord(ifstream& evtfile, unsigned int& variable)
{
    evtfile.read((char*)buffer,BUFFER_BYTES);
    if(evtfile)
    {
        variable = buffer[0];
        return true;
    }

    return false;
}

// read two words from the input file, combine them, and store as a
// a two-word variable
bool readTwoWordVariable(ifstream& evtfile, unsigned int& variable)
{
    unsigned int var1;
    unsigned int var2;

    if(readWord(evtfile, var1) && readWord(evtfile, var2))
    {
        variable = (var2 << 16) | var1;
        return true;
    }

    //cerr << "Error: failed to read a two-word variable from the input file."  << endl;
    return false;
}

bool readEventHeader(ifstream& evtfile)
{
    // start reading the header of the new event
    if(readTwoWordVariable(evtfile, rawEvent.size)
    && readTwoWordVariable(evtfile, rawEvent.evtType)
    && readTwoWordVariable(evtfile, rawEvent.chNo)
    && readTwoWordVariable(evtfile, rawEvent.timetag))
    {
        // successfully read event header
        return true;
    }

    return false;
}

// read the extras data from the input file (only applicable if you're in
// the DPP event body)
bool readExtras(ifstream& evtfile)
{
    readWord(evtfile, rawEvent.extraSelect);
    switch(rawEvent.extraSelect)
    {
        unsigned int dummy;
        case 0:
            readWord(evtfile, dummy); // baseline value * 4 (this information currently discarded)
            rawEvent.baseline = dummy/4;
            readWord(evtfile, rawEvent.extTime);
            return true;

        case 2:
            readWord(evtfile, dummy);
            rawEvent.flags = (dummy & CONFIG_FLAGS_MASK);
            rawEvent.fineTime = (dummy & FINETIME_MASK);
            readWord(evtfile, rawEvent.extTime);
            return true;

        case 5:
            readWord(evtfile, rawEvent.NZC);
            readWord(evtfile, rawEvent.PZC);
            rawEvent.fineTime = ((double)(8192-rawEvent.NZC)/
                (rawEvent.PZC-rawEvent.NZC));
            return true;

        default:
            // other cases not currently implemented
            cerr << "Error: encountered unimplemented value of extraSelect" << endl;
            return false;
    }
}

// read the event's waveform from the input file
bool readWaveformData(ifstream& evtfile)
{
    if(!rawEvent.waveform)
    {
        rawEvent.waveform = new vector<int>;
    }

    readTwoWordVariable(evtfile, rawEvent.nSamp);
    for(int i=0;(size_t)i<rawEvent.nSamp;i++)
    {
        evtfile.read((char*)buffer,BUFFER_BYTES);
        rawEvent.waveform->push_back(buffer[0]);
    }

    return true;
}

// read the entire DPP event body from the input file
bool readDPPEventBody(ifstream& evtfile)
{
    readExtras(evtfile);
    readWord(evtfile, rawEvent.sgQ);
    readWord(evtfile, rawEvent.lgQ);

    // "pile-up rejection" included in datastream, but not implemented in
    // current acquisition software - discard this data
    evtfile.read((char*)buffer,BUFFER_BYTES);

    // "probe" indicates whether an additional diagnostic waveform will be
    // captured along with the normal waveform. The top bit turns on this
    // probe waveform; the bottom two bits describe the probe type.
    // Note that this is a diagnostic tool for acquisition turned off
    // during production runs
    evtfile.read((char*)buffer,BUFFER_BYTES);

    return readWaveformData(evtfile);
}

// read the entire waveform event body from the input file
bool readWaveformEventBody(ifstream& evtfile)
{
    return readWaveformData(evtfile);
}

// read a single event from the input file
bool readEvent(ifstream& evtfile)
{
    // to prepare for reading a new event, reset all DPP-specific event variables
    //rawEvent.extTime = 0;
    //rawEvent.fineTime = 0;
    //rawEvent.sgQ = 0;
    //rawEvent.lgQ = 0;
    if(rawEvent.waveform)
    {
        rawEvent.waveform->clear();
    }

    if(!readEventHeader(evtfile))
    {
        return false;
    }
    
    // start reading event body
    if(rawEvent.evtType==1)
    {
        return readDPPEventBody(evtfile);
    }

    if (rawEvent.evtType==2)
    {
        return readWaveformEventBody(evtfile);
    }

    cerr << "Error: evtType must be either 1 (DPP mode) or 2 (Waveform mode)." << endl;
    return false;
}

void readRawData(string inFileName, string outFileName, vector<string> channelMap)
{
    cout << "Creating " << outFileName << endl;

    // attempt to open input file
    ifstream inFile;
    inFile.open(inFileName,ios::binary);
    if (!inFile)
    {
        cerr << "Failed to open " << inFileName << ". Please check that the file exists" << endl;
        exit(1);
    }
    cout << inFileName << " opened successfully. Start reading events..." << endl;

    // create output file and ROOT tree for storing events
    TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");
    vector<TTree*> orchardSplit; // channel-specific DPP events
    vector<TTree*> orchardSplitW;// channel-specific waveform events

    // Create trees to be filled with sorted data
    // Each channel has a separate tree for DPP data and for waveform mode data
    for(int i=0; (size_t)i<channelMap.size(); i++)
    {
        orchardSplit.push_back(new TTree((channelMap[i]).c_str(),""));
        orchardSplitW.push_back(new TTree((channelMap[i]+"W").c_str(),""));

        if(channelMap[i]=="-")
        {
            continue;
        }

        branchRaw(orchardSplit[i]);
        orchardSplit[i]->SetDirectory(outFile);

        branchRaw(orchardSplitW[i]);
        orchardSplitW[i]->SetDirectory(outFile);
    }

    // create a vector to hold waveform samples
    rawEvent.waveform = new vector<int>;

    int extTime = 0;
    double prevTimetag = 0;
    int prevEvtType = 0;


    // count the number of events processed
    long rawNumberOfEvents = 0;
    long rawNumberOfDPPs = 0;
    long rawNumberOfWaveforms = 0;

    // start looping through the evtfile to extract events
    while(!inFile.eof())
    {
        readEvent(inFile);

        rawNumberOfEvents++;

        // manually increment extended time, if necessary
        /*if((prevTimetag > (double)pow(2,32) - 50000000) &&
            rawEvent.timetag < prevTimetag)
        {
           extTime++; 
        }*/

        // if LED is used instead of CFD, clear the finetime variable
        if(rawEvent.chNo==0 || rawEvent.chNo==2 || rawEvent.chNo==6)
        {
            rawEvent.fineTime = 0;
        }

        /*if(prevEvtType!=rawEvent.evtType)
        {
            extTime=0;
            prevTimetag=0;
        }

        rawEvent.extTime = extTime;
        */
           
        if(rawEvent.evtType==1)
        {
            rawNumberOfDPPs++;
            orchardSplit[rawEvent.chNo]->Fill();
        }

        else if(rawEvent.evtType==2)
        {
            rawNumberOfWaveforms++;
            orchardSplitW[rawEvent.chNo]->Fill();
        }

        prevEvtType = rawEvent.evtType;
        prevTimetag = rawEvent.timetag;

        // print progress every 10000 events
        if (rawNumberOfEvents%10000 == 0)
        {
            cout << "Processed " << rawNumberOfEvents << " events\r";
            fflush(stdout);
        }
    }

    // reached end of input file - print statistics and clean up
    cout << "Finished processing event file" << endl;
    cout << "Total events: " << rawNumberOfEvents << endl;
    cout << "Total number of DPP-mode events processed = " << rawNumberOfDPPs << endl;
    cout << "Total number of waveform-mode events processed = " << rawNumberOfWaveforms << endl;

    inFile.close();
    outFile->Write();
    outFile->Close();
}
