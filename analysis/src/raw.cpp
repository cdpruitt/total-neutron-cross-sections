/******************************************************************************
  raw.cpp
 ******************************************************************************/

// raw.cpp reads raw event data from a binary .evt file. It produces a ROOT file
// with a ROOT tree for each digitizer channel. Each tree contains all the data
// taken on that channel. No other processing of the data is done (basically,
// this method maps the original evt file to several ROOT trees).
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
//   DPP event body (event type 1):
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
//   WAVEFORM event body (event type 2):
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
#include "TH1.h"

#include "../include/dataStructures.h" // defines the C-structs that hold each event's data
#include "../include/branches.h"       // methods to map C-structs that hold raw data to ROOT trees, and vice-versa
#include "../include/softwareCFD.h"    // for improving timing accuracy using waveform data
#include "../include/config.h"

#include "../include/raw.h"            // declarations of functions used for reading raw data

using namespace std;

const unsigned int CONFIG_FLAGS_MASK = 0xfc00;
const unsigned int FINETIME_MASK = 0x03ff;

const unsigned int BUFFER_SIZE = 2; // in bytes
unsigned short buffer[BUFFER_SIZE/(sizeof(unsigned short))]; // for holding words read from the input file

extern Config config;

// read a word from the input file and store in a variable
bool readWord(ifstream& file, unsigned int& variable)
{
    file.read((char*)buffer,BUFFER_SIZE);

    if(file)
    {
        variable = buffer[0];
        return true;
    }

    return false;
}

// discard the next word from the input file
bool discardWord(ifstream& file)
{
    file.read((char*)buffer,BUFFER_SIZE);

    if(file)
    {
        return true;
    }

    return false;
}

// read two words from the input file, combine them, and store as a
// a two-word variable
bool readTwoWords(ifstream& file, unsigned int& variable)
{
    unsigned int var1, var2;

    if(readWord(file, var1) && readWord(file, var2))
    {
        variable = (var2 << 16) | var1;
        return true;
    }

    return false;
}

// read all the components of an event header
bool readEventHeader(ifstream& file, RawEvent& rawEvent)
{
    return (   readTwoWords(file, rawEvent.size)
            && readTwoWords(file, rawEvent.evtType)
            && readTwoWords(file, rawEvent.chNo)
            && readTwoWords(file, rawEvent.timetag));
}

// read the short-gate integrated charge variable
bool readSGQ(ifstream& file, RawEvent& rawEvent)
{
    return readWord(file, rawEvent.sgQ);
}

// read the long-gate integrated charge variable
bool readLGQ(ifstream& file, RawEvent& rawEvent)
{
    return readWord(file, rawEvent.lgQ);
}

// "pile-up rejection" is a bit that indicates whether two events have occurred
// so closely in time that the second one cannot be processed (and should be
// rejected).
// this feature is not yet implemented on the digitizer, so this word is
// discarded
bool readPileUpRejection(ifstream& file, RawEvent& rawEvent)
{
    return discardWord(file);
}

// "probe" indicates whether an additional diagnostic waveform will be
// captured along with the normal waveform. The top bit turns on this
// probe waveform; the bottom two bits describe the probe type.
// this feature is not yet implemented on the digitizer, so this word is
// discarded
bool readProbe(ifstream& file, RawEvent& rawEvent)
{
    return discardWord(file);
}

// read the extras data from the input file
// (only applicable if the event is DPP-mode)
bool readExtras(ifstream& file, RawEvent& rawEvent)
{
    readWord(file, rawEvent.extraSelect);

    switch(rawEvent.extraSelect)
    {
        unsigned int dummy;
        case 0:
            readWord(file, dummy);
            rawEvent.baseline = dummy/4; // The first word is (baseline value) * 4
            readWord(file, rawEvent.extTime);
            return true;

        case 2:
            readWord(file, dummy);
            rawEvent.flags = (dummy & CONFIG_FLAGS_MASK);
            rawEvent.fineTime = (dummy & FINETIME_MASK);
            readWord(file, rawEvent.extTime);
            return true;

        case 5:
            readWord(file, rawEvent.NZC);
            readWord(file, rawEvent.PZC);
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
bool readWaveformData(ifstream& file, RawEvent& rawEvent)
{
    if(readTwoWords(file, rawEvent.nSamp))
    {
        // empty the vector to prepare for filling with new waveform
        rawEvent.waveform.clear();

        for(unsigned int i=0; i<rawEvent.nSamp; i++)
        {
            unsigned int dummy;
            readWord(file, dummy);
            rawEvent.waveform.push_back(dummy);
        }

        return true;
    }

    return false;
}

// read the entire DPP event body from the input file
bool readDPPEventBody(ifstream& file, RawEvent& rawEvent)
{
    return (
          readExtras(file, rawEvent)
       && readSGQ(file, rawEvent)
       && readLGQ(file, rawEvent)
       && readPileUpRejection(file, rawEvent)
       && readProbe(file, rawEvent)
       && readWaveformData(file, rawEvent));
}

// read the entire waveform event body from the input file
bool readWaveformEventBody(ifstream& file, RawEvent& rawEvent)
{
    return readWaveformData(file, rawEvent);
}

// read a single event from the input file
bool readEvent(ifstream& file, RawEvent& rawEvent)
{
    if(!readEventHeader(file, rawEvent))
    {
        return false;
    }
    
    if(rawEvent.evtType==1)
    {
        // DPP mode
        return readDPPEventBody(file, rawEvent);
    }

    if(rawEvent.evtType==2)
    {
        // waveform mode
        return readWaveformEventBody(file, rawEvent);
    }

    cerr << "Error: evtType must be either 1 (DPP mode) or 2 (Waveform mode)." << endl;
    return false;
}

int readRawData(string inFileName, string outFileName, string DPPTreeName, string WaveformTreeName)
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

    TTree* DPPTree = new TTree(DPPTreeName.c_str(),"");
    TTree* WaveformTree = new TTree(WaveformTreeName.c_str(),"");

    RawEvent rawEvent; // for holding raw event data from the input file in preparation for transfer to a ROOT tree

    DPPTree->Branch("completeTime",&rawEvent.completeTime,"completeTime/d");
    DPPTree->Branch("fineTime",&rawEvent.fineTime,"fineTime/d");
    DPPTree->Branch("cycleNumber",&rawEvent.cycleNumber,"cycleNumber/i");
    DPPTree->Branch("chNo",&rawEvent.chNo,"chNo/i");
    DPPTree->Branch("extTime",&rawEvent.extTime,"extTime/i");
    DPPTree->Branch("timetag",&rawEvent.timetag,"timetag/i");
    DPPTree->Branch("sgQ",&rawEvent.sgQ,"sgQ/i");
    DPPTree->Branch("lgQ",&rawEvent.lgQ,"lgQ/i");
    DPPTree->Branch("baseline",&rawEvent.baseline,"baseline/i");
    DPPTree->Branch("waveform",&rawEvent.waveform);

    WaveformTree->Branch("completeTime",&rawEvent.completeTime,"completeTime/d");
    WaveformTree->Branch("fineTime",&rawEvent.fineTime,"fineTime/d");
    WaveformTree->Branch("cycleNumber",&rawEvent.cycleNumber,"cycleNumber/i");
    WaveformTree->Branch("chNo",&rawEvent.chNo,"chNo/i");
    WaveformTree->Branch("extTime",&rawEvent.extTime,"extTime/i");
    WaveformTree->Branch("timetag",&rawEvent.timetag,"timetag/i");
    WaveformTree->Branch("sgQ",&rawEvent.sgQ,"sgQ/i");
    WaveformTree->Branch("lgQ",&rawEvent.lgQ,"lgQ/i");
    WaveformTree->Branch("baseline",&rawEvent.baseline,"baseline/i");
    WaveformTree->Branch("waveform",&rawEvent.waveform);

    // count the number of events processed
    long rawNumberOfEvents = 0;
    long rawNumberOfDPPs = 0;
    long rawNumberOfWaveforms = 0;

    unsigned int prevEvtType = 0;
    rawEvent.cycleNumber = 0;

 
    // start looping through the evtfile to extract events
    while(!inFile.eof())
    {
        if(!readEvent(inFile, rawEvent))
        {
            cerr << "Failed to read event number " << rawNumberOfEvents << ". Ending read of raw events." << endl;

            inFile.close();
            outFile->Write();
            outFile->Close();

            return 1;
        }

        // Assign a time to each event
        rawEvent.completeTime =
            double(pow(2,31)*rawEvent.extTime) +
            rawEvent.timetag; // in samples
        rawEvent.completeTime *= config.digitizer.SAMPLE_PERIOD; // converts from samples to ns

        // Discard events with digitizer rollover error (i.e., incrementing extTime before clearing the timetag, near the rollover period at 2^31 bits)
/*
        if(macropulseEvent.timeExtTime > prevExtTime
                && macropulseEvent.timeTimetag > pow(2,31)-1000)
        {
            cerr << "Error: digitizer rollover error. Skipping event. macroNo = " << macropulseEvent.macroNo << ", extTime = " << macropulseEvent.timeExtTime << ", timetag = " << macropulseEvent.timeTimetag << endl;
            cerr << "Skipping to next target changer event..." << endl;

            i++;
            continue;
        }
        */

        // sync times to the macropulse timing channel by applying an channel-dependent offset (accounts for cable delay)
        switch(rawEvent.chNo)
        {
            case 0:
                rawEvent.completeTime += config.timeOffsets.TARGET_CHANGER_TIME_OFFSET;
                break;
            case 1:
                break;
            case 2:
                rawEvent.completeTime += config.timeOffsets.MONITOR_TIME_OFFSET;
                break;
            case 4:
                rawEvent.completeTime += config.timeOffsets.DETECTOR_TIME_OFFSET;

                // use CFD to improve timing precision
                rawEvent.fineTime = calculateCFDTime(
                            rawEvent.waveform,
                            rawEvent.baseline,
                            config.softwareCFD.CFD_FRACTION,
                            config.softwareCFD.CFD_DELAY); // CFD time in samples
                if(rawEvent.fineTime>=0)
                {
                    // recovered a good fine time for this event
                    rawEvent.completeTime += (rawEvent.fineTime-config.softwareCFD.CFD_TIME_OFFSET)*config.digitizer.SAMPLE_PERIOD;
                }
                    break;
            case 5:
                rawEvent.completeTime += config.timeOffsets.VETO_TIME_OFFSET;
                break;
            case 6:
                rawEvent.completeTime += config.timeOffsets.HIGH_T_DET_TIME_OFFSET;

                // use CFD to improve timing precision
                rawEvent.fineTime = calculateCFDTime(
                            rawEvent.waveform,
                            rawEvent.baseline,
                            config.softwareCFD.CFD_FRACTION,
                            config.softwareCFD.CFD_DELAY); // CFD time in samples
                if(rawEvent.fineTime>=0)
                {
                    // recovered a good fine time for this event
                    rawEvent.completeTime += (rawEvent.fineTime-config.softwareCFD.CFD_TIME_OFFSET)*config.digitizer.SAMPLE_PERIOD;
                }
                    break;
            case 7:
            default:
                cerr << "Error: encountered unimplemented channel number " << rawEvent.chNo << " during complete time assignment. Ending raw data read-in..." << endl;
                return 1;
        }

        if(rawEvent.evtType==1)
        {
            if(prevEvtType==2)
            {
                // first event after a waveform->DPP mode change; increment
                // cycle counter
                rawEvent.cycleNumber++;
            }

            DPPTree->Fill();
            rawNumberOfDPPs++;
        }

        else if(rawEvent.evtType==2)
        {
            WaveformTree->Fill();
            rawNumberOfWaveforms++;
        }

        else
        {
            cerr << "Error: event " << rawNumberOfEvents << "had event type other than 1 (DPP mode) or 2 (waveform mode). Exiting..." << endl;

            inFile.close();
            outFile->Write();
            outFile->Close();

            return(1);
        }

        // print progress every 10000 events
        if (rawNumberOfEvents%10000 == 0)
        {
            cout << "Processed " << rawNumberOfEvents << " events\r";
            fflush(stdout);
        }

        rawNumberOfEvents++;
        prevEvtType = rawEvent.evtType;
    }

    // reached end of input file - print statistics and clean up
    cout << endl;
    cout << "Finished processing event file." << endl;
    cout << "Total events: " << rawNumberOfEvents << endl;
    cout << "Total number of DPP-mode events processed = " << rawNumberOfDPPs << endl;
    cout << "Total number of waveform-mode events processed = " << rawNumberOfWaveforms << endl;

    inFile.close();
    outFile->Write();
    outFile->Close();

    return 0;
}
