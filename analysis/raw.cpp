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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"

using namespace std;

const double SAMPLE_PERIOD = 2; // time of one sample, in ns

// To populate events from the raw file into a ROOT tree, we provide
// an event structure listing all the variables that comprise an event.
// When a new event is added to the tree, a subset of these variables is
// included in the new event being saved.
struct Event
{
    // Event header variables:
    unsigned int size;    // size of event, in bytes
    unsigned int evtType; // "event type", either 1 (DPP) or 2 (waveform)
    unsigned int chNo;    // "channel number" indicates event origin (i.e.,
    // detector, monitor, target changer)
    double timetag;       // "coarse timestamp of event" includes 32 bits of time
    // information. Units are the same as the sample period of
    // the digitizer (e.g. 2 ns or 5 ns)

    // Event body variables:
    unsigned int extraSelect; // indicates the meaning of the "extras" words
    unsigned int sgQ;     // "short gate integrated charge" provides the charge
    // integral over an adjustable range of the event's peak 
    unsigned int lgQ;     // "long gate integrated charge" provides the charge
    // integral over an adjustable range of the event's peak 
    unsigned int nSamp; // number of samples in the event's waveform
    vector<int> waveform; // "digital waveform of event" is a series of waveform
    // samples for each event

    // Variables extracted from "extras"

    unsigned int baseline;
    unsigned int flags;
    unsigned int extTime; // "extended timestamp of event" extends timetag with
    // 16 additional bits for times greater than 2^32
    // sample periods.
    unsigned int fineTime;// "fine timestamp of event" sub-divides timetag with
    // 10 additional bits of time granularity, with units of
    // (sample period)/2^10 units.
} ev;

// Raw data is stored as hexadecimal words (16 bits long) in the .evt files
int const BufferWords = 1; // number of chars per buffer word
int const BufferBytes = BufferWords*2; // number of bytes per buffer word

// track number of processed events of each type
struct Statistics
{
    long numberOfEvents = 0;
    long numberOfDPPs = 0;
    long numberOfWaveforms = 0;
    long numberOfCh0Waveforms = 0;
    long numberOfCh2Waveforms = 0;
    long numberOfCh4Waveforms = 0;
} stats;

// extract a single event from the raw event file and fill its data into the
// tree
void readEvent(ifstream& evtfile)
{
    // clear all DPP-specific event variables to prepare for reading event
    ev.extTime = 0;
    ev.fineTime = 0;
    ev.sgQ = 0;
    ev.lgQ = 0;

    unsigned short buffer[BufferWords];

    // start reading event header (common to all events)

    evtfile.read((char*)buffer,BufferBytes);
    unsigned short size1 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    unsigned short size2 = buffer[0];
    ev.size = (size2 << 16) | size1;

    evtfile.read((char*)buffer,BufferBytes);
    unsigned short evtType1 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    unsigned short evtType2 = buffer[0];
    ev.evtType = (evtType2 << 16) | evtType1;

    evtfile.read((char*)buffer,BufferBytes);
    unsigned short chNo1 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    unsigned short chNo2 = buffer[0];
    ev.chNo = (chNo2 << 16) | chNo1;

    evtfile.read((char*)buffer,BufferBytes);
    unsigned short timetag1 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    unsigned short timetag2 = buffer[0];
    ev.timetag = (timetag2<< 16) | timetag1;
    ev.timetag *= SAMPLE_PERIOD; // timetag converted from samples to ns
    // finished reading event header

    if(ev.evtType==1)
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

        evtfile.read((char*)buffer,BufferBytes);
        ev.extraSelect = buffer[0];

        evtfile.read((char*)buffer,BufferBytes);
        unsigned int extras1 = buffer[0];

        evtfile.read((char*)buffer,BufferBytes);
        unsigned int extras2 = buffer[0];

        switch(ev.extraSelect)
        {
            case 2:
                // retrieve extended time from bits 16-31 (extras2)
                ev.extTime = extras2;
                // retrieve configuration file flags from bits 10:15 (0xfc00)
                ev.flags = (extras1 & 0xfc00);
                // retrieve fine time from bits 0:9 (0x03ff)
                ev.fineTime = (extras1 & 0x03ff);
                break;

                // other cases not currently implemented
            default:
                cout << "Error: extraSelect != 2; other values of extraSelect not implemented" << endl;
                exit(1);
                break;
        }

        evtfile.read((char*)buffer,BufferBytes);
        ev.sgQ = buffer[0];

        evtfile.read((char*)buffer,BufferBytes);
        ev.lgQ = buffer[0];

        // "pile-up rejection" flag (not implemented in current acquisition software,
        // so discard this data)
        //puRej = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);

        // "probe" indicates whether an additional diagnostic waveform will be captured along with the
        // normal waveform. The top bit turns on this probe waveform; the bottom two bits describe the
        // probe type.
        // Note that this is a diagnostic tool for acquisition turned off during production runs
        //probe = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);

        stats.numberOfDPPs++;
    }

    else if (ev.evtType==2)
    {
        // waveform mode event
        stats.numberOfWaveforms++;
    }

    else
    {
        cout << "Error: evtType must be either 1 (DPP mode) or 2 (Waveform mode)." << endl;
        exit(1);
    }

    // the DPP-mode only data has been read out, so now let's read out the data
    // that's always present (the number of waveform samples and the waveforms)

    evtfile.read((char*)buffer,BufferBytes);
    unsigned short nSamp1 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    unsigned short nSamp2 = buffer[0];
    ev.nSamp = (nSamp2 << 16) | nSamp1;

    ev.waveform.clear();
    for(int i=0;(size_t)i<ev.nSamp;i++)
    {
        evtfile.read((char*)buffer,BufferBytes);
        ev.waveform.push_back(buffer[0]);
    }


    // Update statistics on events

    stats.numberOfEvents++;
    if(ev.evtType==2)
    {
        if(ev.chNo==0)
        {
            stats.numberOfCh0Waveforms++;
        }

        if(ev.chNo==2)
        {
            stats.numberOfCh2Waveforms++;
        }

        if(ev.chNo==4)
        {
            stats.numberOfCh4Waveforms++;
        }
    }
}

int main(int argc, char* argv[])
{
    cout << endl << "Entering raw sort..." << endl;

    string inFileName = argv[1];
    string outFileName = argv[2];

    // Open a ROOT file to store the event data
    TFile *outFile;
    outFile = new TFile(outFileName.c_str(),"UPDATE");

    if(outFile->Get("tree"))
    {
        cout << "Found previously existing raw sort " << outFileName << ". Skipping raw sort." << endl;
        exit(0);
    }

    // Create a ROOT tree for storing events 
    TTree* tree = new TTree("tree","");
    tree->Branch("evtType",&ev.evtType,"evtType/i");
    tree->Branch("chNo",&ev.chNo,"chNo/i");
    tree->Branch("extTime",&ev.extTime,"extTime/i");
    tree->Branch("timetag",&ev.timetag,"timetag/d");
    tree->Branch("fineTime",&ev.fineTime,"fineTime/i");
    tree->Branch("sgQ",&ev.sgQ,"sgQ/i");
    tree->Branch("lgQ",&ev.lgQ,"lgQ/i");
    tree->Branch("waveform",&ev.waveform);

    // attempt to open input file
    ifstream inFile;
    inFile.open(inFileName,ios::binary);
    if (!inFile)
    {
        cout << "Failed to open " << inFileName << ". Please check that the file exists" << endl;
        exit(1);
    }

    cout << inFileName << " opened successfully. Start reading events..." << endl;

    // we're now pointing at the first 16-bit word in the data stream
    // start looping through the evtfile to extract events
    while(!inFile.eof() /* use to truncate sort && stats.numberOfEvents<1000000*/)
    {
        readEvent(inFile);

        // Add event to tree
        tree->Fill();

        if (stats.numberOfEvents%10000 == 0)
        {
            cout << "Processed " << stats.numberOfEvents << " events\r";
            fflush(stdout);
        }
    }

    // reached end of input file
    cout << "Finished processing event file" << endl;
    cout << "Total events: " << stats.numberOfEvents << endl;
    cout << "Total number of DPP-mode events processed = " << stats.numberOfDPPs << endl;
    cout << "Total number of waveform-mode events processed = " << stats.numberOfWaveforms << endl;

    inFile.close();

    outFile->Write();
    outFile->Close();
}
