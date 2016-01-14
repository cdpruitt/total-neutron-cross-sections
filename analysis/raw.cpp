// Reads output from washudaq-X.X and sorts into ROOT and text files
// 
// The raw event file structure should be as follows:
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
//      Extra select    |   uint16; describes the contents of Extras
//      Extras          |   uint32; additional PSD data (see documentation)
//      Zero Crossing   |   uint16; interpolated zero-crossing, in picoseconds
//      Short Gate Q    |   uint16; charge integrated in the short gate
//      Long Gate Q     |   uint16; charge integrated in the long gate
//      Baseline        |   uint16; baseline level of ADC
//      Pile up rej.    |   uint16; pile-up rejection flag (not yet implemented)
//      Probe info      |   uint16; turns on an analog probe for examining waveform
//      Num. of samples |   uint32; number of waveform samples collected in Samples
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
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"

using namespace std;

const string analysispath =  "/media/Drive3/";

struct event
{
    unsigned int chNo; // describe data stream origin (i.e., detector, monitor, target changer)
    unsigned int evtType; // either 1 (DPP) or 2 (waveform)

    double timetag; // 1 sample granularity, 32 bits
    unsigned int extTime, fineTime; // provide additional bits of granularity 

    unsigned int sgQ, lgQ; // integrated charge, from a short gate or a long gate

    vector<int> waveform; // contains all waveform samples for each event to allow for corrections in analysis
} ev;

// Buffer variables
int const BufferWords = 1;
int const BufferBytes = BufferWords*2;
unsigned short buffer[BufferWords];
unsigned short *point;

// data variables to read from input file 
unsigned int size, evtType, chNo, sgQ, lgQ, extTime, fineTime, nSamp, extraSelect, extras1, extras2;
double timetag;
vector<int> waveform; // holds waveform samples for each event

// Tree definition
TTree* tree;

void readEvent(ifstream& evtfile)
{
    // start reading header

    // size is the number of bytes in the event (self-inclusive)
    unsigned short size1 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    unsigned short size2 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    size = (size2 << 16) | size1;

    // evtType is either 1 (DPP mode) or 2 (waveform mode)
    unsigned short evtType1 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    unsigned short evtType2 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    evtType = (evtType2 << 16) | evtType1;

    // chNo ranges from 0-7
    unsigned short chNo1 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    unsigned short chNo2 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    chNo = (chNo2 << 16) | chNo1;

    // timetag is the coarse trigger time for this event, in 2ns increments
    unsigned short timetag1 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    unsigned short timetag2 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    timetag = (timetag2<< 16) | timetag1;
    timetag *= 2; // timetag converted from samples to ns

    if(evtType==1)
    {
        // DPP EVENT BODY unpacking 

        // EXTRAS is a 32-bit word whose content varies with the value of EXTRA_SELECT.
        // [DEFAULT]    0: extended timestamp (bits 16-31) and baseline*4 (bits 0-15)
        //              1: extended timestamp (bits 16-31) and flags (bits 0-15?)
        //              2: extended timestamp (bits 16-31), flags (bits 10-15), and fine timestamp
        //                  (bits 0-9)
        //              3: pulse peak value (bits 0-15) [UNTESTED FEATURE]
        //              5: CFD positive ZC (bits 16-31) and negative ZC (bits 0-15)
        //              7: fixed value of 0x12345678

        extraSelect = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);

        extras1 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);

        extras2 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);

        switch(extraSelect)
        {
            case 2:
                // retrieve extended time from bits 16-31 (extras2)
                extTime = extras2;
                // retrieve fine time from bits 0:9 (0x03ff)
                fineTime = (extras1 & 0x03ff);
                break;

        }

        // sgQ is the short gate integrated charge, in digitizer units 
        sgQ = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);

        // lgQ is the long gate integrated charge, in digitizer units 
        lgQ = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);

        // puRej is a pile-up rejection flag (not used in current acquisition
        // software - discard this data)
        //puRej = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);

        // probe indicates whether an additional analog waveform will be captured along with the
        // input trace. The top bit turns on the analog probe; the bottom two bits describe the
        // probe type.
        //probe = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
    }

    // the DPP-mode only data has been read out, so now let's read out the data
    // that's present in both DPP and waveform modes (number of samples and the
    // waveform data)

    // nSamp is the number of waveform samples that follow (in LIST mode, nSamp = 0)
    unsigned short nSamp1 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    unsigned short nSamp2 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    nSamp = (nSamp2 << 16) | nSamp1;

    // read waveform, consisting of nSamp samples
    waveform.clear();
    
    for(int i=0;i<nSamp;i++)
    {
        waveform.push_back(buffer[0]);
        evtfile.read((char*)buffer,BufferBytes);
    }

    // fill tree w/ event data
    ev.chNo = chNo;
    ev.evtType = evtType;
    ev.extTime = extTime;
    ev.timetag = timetag;
    ev.fineTime = fineTime;
    ev.sgQ = sgQ;
    ev.lgQ = lgQ;
    ev.waveform = waveform; 

    // ignore 'junk' macropulse signals with off-scale lgQ
    // --> not sure what these events are... didn't notice them during the
    // experiment
    if (lgQ != 65535 || chNo !=0)
    {
        tree->Fill();
    }

    // clear all DPP-specific event variables to prepare for next event
    extTime = 0;
    fineTime = 0;
    sgQ = 0;
    lgQ = 0;
}

void processRun(string evtname)
{
    // attempt to process the event file 'evtname' as defined in main()
    ifstream evtfile;
    evtfile.open(evtname,ios::binary);

    if (!evtfile)
    {
        cout << "Failed to open " << evtname << ". Please check that the file exists" << endl;
        abort();
    }

    else // successfully opened event file; start processing events
    {
        cout << evtname << " opened successfully. Start reading events..." << endl;

        evtfile.read((char*)buffer,BufferBytes);

        point = buffer;
        // we're now pointing at the first 16-bit word in the data stream

        // start looping through the evtfile to extract events
        int nE = 0;

        while(!evtfile.eof() /* use to truncate sort */ && nE<6000000)
        {

            readEvent(evtfile); // extract raw data from event file

            nE++;

            if (nE%10000 == 0)
            {
                cout << "Processed " << nE << " events\r";
                fflush(stdout);
            }
        }

        // Input file finished
        cout << "Finished processing event file" << endl;
        cout << "Total events: " << nE << endl;
    }

    evtfile.close();
}

int main(int argc, char* argv[])
{
    // read in data run location
    string runDir = argv[1];
    string runNo = argv[2];

    // Create a tree for this run
    TFile *file;

    stringstream treeName;
    stringstream fileName;
    treeName << runDir << "-" << runNo; 
    fileName << analysispath <<"analysis/" << runDir << "/" << treeName.str() << ".root";

    file = new TFile(fileName.str().c_str(),"UPDATE");

    if(file->Get("tree"))
    {
        tree = (TTree*)file->Get("tree");
        cout << "Found previously existing tree - skipping sorting " << treeName << endl;
    }

    else
    {
        tree = new TTree("tree","");
        cout << "Created ROOT tree " << treeName.str() << endl;

        tree->Branch("evtType",&ev.evtType,"evtType/i");
        tree->Branch("chNo",&ev.chNo,"chNo/i");
        tree->Branch("extTime",&ev.extTime,"extTime/i");
        tree->Branch("timetag",&ev.timetag,"timetag/d");
        tree->Branch("fineTime",&ev.fineTime,"fineTime/i");
        tree->Branch("sgQ",&ev.sgQ,"sgQ/i");
        tree->Branch("lgQ",&ev.lgQ,"lgQ/i");
        tree->Branch("waveform",&ev.waveform);
    }

    stringstream runName;
    runName << analysispath <<"output/" << runDir << "/data-" << runNo << ".evt";
    processRun(runName.str());

    file->Write();

    file->Close();

}
