#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "../include/raw.h"
#include "../include/dataStructures.h"
#include "../include/physicalConstants.h"

extern RawEvent rawEvent;

/******************************************************************************
  raw.cpp
 ******************************************************************************/

// This is called by ./sort.sh and is the first step in analyzing the total
// cross section data.
// It reads raw hexidecimal data from an .evt file and creates a ROOT file,
// runX-Y_raw.root (where X is the run number, and Y is the sub-run number),
// with one ROOT tree containing all events from the input .evt file.
// 
// This file expects .evt file data with the following structure:
//
//      ---------------------------------------------------------------------------
//      EVENT HEADER:
//
//      Size            |   uint32; number of bytes in the event (self-inclusive)
//      Event Type      |   uint32; possible values are 1 (DPP) or 2 (waveform)
//      Channel         |   uint32; channel number of event (range: 0-7)
//      Time Tag        |   uint32; coarse time (2 ns units) of event trigger
//
//      (in every event, there is an EVENT HEADER)
//      ---------------------------------------------------------------------------
//
//      ---------------------------------------------------------------------------
//      DPP EVENT BODY:
//
//      Extra select    |   uint16; describes the contents of Extras
//      Extras          |   uint32; additional PSD data (see DPP Event body section)
//      Short Gate Q    |   uint16; charge integrated in the short gate
//      Long Gate Q     |   uint16; charge integrated in the long gate
//      Baseline        |   uint16; baseline level of ADC
//      Pile up rej.    |   uint16; pile-up rejection flag (not yet implemented)
//      Probe info      |   uint16; turns on an analog probe for examining waveform
//      Num. of samples |   uint32; number of waveform samples collected in Samples
//      Samples         |   series of uint16 giving raw waveform samples
//
//      (only events of type 1 have DPP EVENT BODY)
//      ---------------------------------------------------------------------------
//
//      ---------------------------------------------------------------------------
//
//      WAVEFORM EVENT BODY
//
//      Num. of samples |   uint32; number of waveform samples to follow (mixed mode)
//      Samples         |   series of uint16 giving raw waveform samples
//
//      (only events of type 2 have WAVEFORM EVENT BODY)
//      ---------------------------------------------------------------------------

using namespace std;

// extract a single event from the raw event file and fill its data into the
// tree
int readEvent(ifstream& evtfile)
{
    // clear all DPP-specific event variables to prepare for reading event
    rawEvent.extTime = 0;
    rawEvent.fineTime = 0;
    rawEvent.sgQ = 0;
    rawEvent.lgQ = 0;

    unsigned short buffer[BUFFER_WORDS];
    // we're now pointing at the first 16-bit word in the data stream

    // start reading event header (common to all events)

    evtfile.read((char*)buffer,BUFFER_BYTES);
    unsigned short size1 = buffer[0];
    evtfile.read((char*)buffer,BUFFER_BYTES);
    unsigned short size2 = buffer[0];
    rawEvent.size = (size2 << 16) | size1;

    evtfile.read((char*)buffer,BUFFER_BYTES);
    unsigned short evtType1 = buffer[0];
    evtfile.read((char*)buffer,BUFFER_BYTES);
    unsigned short evtType2 = buffer[0];
    rawEvent.evtType = (evtType2 << 16) | evtType1;

    evtfile.read((char*)buffer,BUFFER_BYTES);
    unsigned short chNo1 = buffer[0];
    evtfile.read((char*)buffer,BUFFER_BYTES);
    unsigned short chNo2 = buffer[0];
    rawEvent.chNo = (chNo2 << 16) | chNo1;

    evtfile.read((char*)buffer,BUFFER_BYTES);
    unsigned short timetag1 = buffer[0];
    evtfile.read((char*)buffer,BUFFER_BYTES);
    unsigned short timetag2 = buffer[0];
    rawEvent.timetag = (timetag2<< 16) | timetag1;
    rawEvent.timetag *= SAMPLE_PERIOD; // timetag converted from samples to ns
    // finished reading event header

    if(rawEvent.evtType==1)
    {
        // start reading event data specific to DPP events

        // EXTRAS is a 32-bit word whose content varies with the value of EXTRA_SELECT.
        // [DEFAULT]    0: extended timestamp (bits 16-31) and baseline*4 (bits 0-15)
        //              1: extended timestamp (bits 16-31) and flags (bits 0-15?)
        //              2: extended timestamp (bits 16-31), flags (bits 10-15), and fine timestamp
        //                  (bits 0-9)
        //              3: pulse peak value (bits 0-15) [UNTESTED FEATURE]
        //              5: CFD positive ZC (bits 16-31) and negative ZC (bits 0-15)
        //              7: fixed value of 0x12345678

        evtfile.read((char*)buffer,BUFFER_BYTES);
        rawEvent.extraSelect = buffer[0];

        evtfile.read((char*)buffer,BUFFER_BYTES);
        unsigned int extras1 = buffer[0];

        evtfile.read((char*)buffer,BUFFER_BYTES);
        unsigned int extras2 = buffer[0];

        switch(rawEvent.extraSelect)
        {
            case 2:
                // retrieve extended time from bits 16-31 (extras2)
                rawEvent.extTime = extras2;
                // retrieve configuration file flags from bits 10:15 (0xfc00)
                rawEvent.flags = (extras1 & 0xfc00);
                // retrieve fine time from bits 0:9 (0x03ff)
                rawEvent.fineTime = (extras1 & 0x03ff);
                break;

                case 5:
                // retrieve positive and negative zero-crossings
                rawEvent.PZC = extras2;
                rawEvent.NZC = extras1;

                // other cases not currently implemented
            default:
                cout << "Error: encountered unimplemented value of extraSelect" << endl;
                return 1;
        }

        evtfile.read((char*)buffer,BUFFER_BYTES);
        rawEvent.sgQ = buffer[0];

        evtfile.read((char*)buffer,BUFFER_BYTES);
        rawEvent.lgQ = buffer[0];

        // "pile-up rejection" flag (not implemented in current acquisition software,
        // so discard this data)
        //puRej = buffer[0];
        evtfile.read((char*)buffer,BUFFER_BYTES);

        // "probe" indicates whether an additional diagnostic waveform will be captured along with the
        // normal waveform. The top bit turns on this probe waveform; the bottom two bits describe the
        // probe type.
        // Note that this is a diagnostic tool for acquisition turned off during production runs
        //probe = buffer[0];
        evtfile.read((char*)buffer,BUFFER_BYTES);

        rawNumberOfDPPs++;
    }

    else if (rawEvent.evtType==2)
    {
        // waveform mode event
        rawNumberOfWaveforms++;
    }

    else
    {
        cout << "Error: evtType must be either 1 (DPP mode) or 2 (Waveform mode)." << endl;
        return 1;
    }

    // the DPP-mode only data has been read out, so now let's read out the data
    // that's always present (the number of waveform samples and the waveforms)

    evtfile.read((char*)buffer,BUFFER_BYTES);
    unsigned short nSamp1 = buffer[0];
    evtfile.read((char*)buffer,BUFFER_BYTES);
    unsigned short nSamp2 = buffer[0];
    rawEvent.nSamp = (nSamp2 << 16) | nSamp1;

    rawEvent.waveform->clear();
    for(int i=0;(size_t)i<rawEvent.nSamp;i++)
    {
        evtfile.read((char*)buffer,BUFFER_BYTES);
        rawEvent.waveform->push_back(buffer[0]);
    }


    // Update statistics on events

    rawNumberOfEvents++;
    if(rawEvent.evtType==2)
    {
        if(rawEvent.chNo==0)
        {
            rawNumberOfCh0Waveforms++;
        }

        if(rawEvent.chNo==2)
        {
            rawNumberOfCh2Waveforms++;
        }

        if(rawEvent.chNo==4)
        {
            rawNumberOfCh4Waveforms++;
        }
    }
    return 0;
}

void extractRawData(string inFileName, string outFileName)
{
    cout << endl << "Entering raw sort..." << endl;

    // attempt to open input file
    ifstream inFile;
    inFile.open(inFileName,ios::binary);
    if (!inFile)
    {
        cout << "Failed to open " << inFileName << ". Please check that the file exists" << endl;
        exit(1);
    }

    cout << inFileName << " opened successfully. Start reading events..." << endl;

    TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

    //RawEvent* rawEvent = RawEvent::getInstance();

    // Create a ROOT tree for storing events 
    TTree* tree = new TTree("tree","");
    tree->Branch("evtType",&rawEvent.evtType,"evtType/i");
    tree->Branch("chNo",&rawEvent.chNo,"chNo/i");
    tree->Branch("extTime",&rawEvent.extTime,"extTime/i");
    tree->Branch("timetag",&rawEvent.timetag,"timetag/d");
    tree->Branch("fineTime",&rawEvent.fineTime,"fineTime/i");
    tree->Branch("sgQ",&rawEvent.sgQ,"sgQ/i");
    tree->Branch("lgQ",&rawEvent.lgQ,"lgQ/i");
    tree->Branch("waveform",&rawEvent.waveform);

    // start looping through the evtfile to extract events
    while(!inFile.eof() /* use to truncate sort && rawNumberOfEvents<1000000*/)
    {
        if(readEvent(inFile))
        {
            // readEvent returned an error; end sort
            break;
        }

        // Add event to tree
        tree->Fill();

        if (rawNumberOfEvents%10000 == 0)
        {
            cout << "Processed " << rawNumberOfEvents << " events\r";
            fflush(stdout);
        }
    }

    // reached end of input file
    cout << "Finished processing event file" << endl;
    cout << "Total events: " << rawNumberOfEvents << endl;
    cout << "Total number of DPP-mode events processed = " << rawNumberOfDPPs << endl;
    cout << "Total number of waveform-mode events processed = " << rawNumberOfWaveforms << endl;

    inFile.close();

    outFile->Write();
    outFile->Close();
}
