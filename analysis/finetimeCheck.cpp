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
#include <sstream>
#include "TH1S.h"
#include "TH2S.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectoryFile.h"

using namespace std;

// Buffer variables
int const BufferWords = 1;
int const BufferBytes = BufferWords*2;
unsigned short buffer[BufferWords];
unsigned short *point;

// Header fields for events
unsigned long size;
unsigned long evtype;
unsigned long channel;
unsigned long timetag, timetagP, timeCoin;
unsigned int channelCoin = 1; // used to keep track of coincidence between channels

float finetime = 0;
// the fine time is provided when EXTRA_SELECT=2. finetimeL/finetimeR are assigned to compare
// the L/R dependence of timing.

unsigned int nE = 0; // counter for the total number of events
unsigned int nWavelets = 0; // counter for the number of short MIXED mode waveforms
unsigned int nCFDs = 0; // counter for the number of MIXED mode CFD traces 
unsigned int nBaselines = 0; // counter for the number of MIXED mode Baseline traces
unsigned int nWaveforms = 0; // counter for the number of long MIXED mode waveforms

// Create a ROOT tree for holding DPP data 
TTree* tree;

struct treeEvent {
    unsigned int timetag; // timetag labels the event, and numCh is the number of channels that triggered at that time
    float finetimeL, finetimeR; // contains fine times for detL and detR
    unsigned int lgQL, lgQR; // for doing an energy selection on events
} te;

static unsigned int const COINCIDENCE_REQ = 2;

// Create histograms for DPP data
TH1S* outFT;
TH1S* outCT;
TH2S* outFTLR;

// read EVENT HEADER to determine where the event data should go (using channel number)
// 'out' points to a channel-specific text file
int readHeader(ifstream& evtfile)
{
    // start reading header

    // size is the number of bytes in the event (self-inclusive)
    unsigned short size1 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    unsigned short size2 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    size = (size2 << 16) | size1;

    // evtype is either 1 (DPP mode) or 2 (waveform mode)
    unsigned short evtype1 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    unsigned short evtype2 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    evtype = (evtype2 << 16) | evtype1;

    // channel ranges from 0-7; each channel is considered an independent event generator
    unsigned short channel1 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    unsigned short channel2 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    channel = (channel2 << 16) | channel1;

    // timetag is the coarse trigger time for this event, in 2ns increments
    unsigned short timetag1 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    unsigned short timetag2 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    timetag = (timetag2<< 16) | timetag1;

    return channel;

}


// unpacks EVENT BODY data into ROOT histograms and text output
// 'out' points to a channel-specific text file
void unpack(ifstream& evtfile)
{
    ostringstream histName;
        
    if(evtype==1)
    {
        // DPP EVENT BODY unpacking 

        // EXTRAS is a 32-bit word whose content varies with the value of EXTRA_SELECT.
        // [DEFAULT]    0: extended timestamp (bits 16-31) and baseline*4 (bits 0-15)
        //              1: extended timestamp (bits 16-31) and flags (bits 0-15?)
        //
        //              2: extended timestamp (bits 16-31), flags (bits 10-15), and fine timestamp
        //                  (bits 0-9)
        //              3: pulse peak value (bits 0-15) [UNTESTED FEATURE]
        //              5: CFD positive ZC (bits 16-31) and negative ZC (bits 0-15)
        //              7: fixed value of 0x12345678

        unsigned short extraSelect = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);

        unsigned short extras1 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);

        unsigned short extras2 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        //const unsigned int extras = (extras2 << 16) | extras1;

        stringstream temp;
        if (extraSelect==2)
        {
            // fine time from bits 0:9 (0x03ff)
            finetime = (extras1 & 0x03ff);
            outFT->Fill(finetime);
        }

        else
        {
            cout << "EXTRA_SELECT is != 2, so no fine time available. Please select a run with EXTRA_SELECT = 2" << endl;
        }

        // sgQ is the short gate integrated charge, in digitizer units 
        unsigned short sgQ = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);

        // lgQ is the long gate integrated charge, in digitizer units 
        unsigned short lgQ = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);

        // puRej is a pile-up rejection flag (1 if pile-up detected? 0 if not?)
        unsigned short puRej = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);

        // probe indicates whether an additional analog waveform will be captured along with the
        // input trace. The top bit turns on the analog probe; the bottom two bits describe the
        // probe type.
        unsigned short probe = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        
        // nSamp is the number of waveform samples that follow (in LIST mode, this is 0)
        unsigned short nSamp1 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        unsigned short nSamp2 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        const unsigned int nSamp = (nSamp2 << 16) | nSamp1;

        if(nSamp > 0)
            // We must be in MIXED mode; time to output wavelet.
        {

            for(int i=0;i<nSamp;i++)
            {
                    
                evtfile.read((char*)buffer,BufferBytes);

            }
            
            // done with wavelet; increment the wavelet counter to get ready for the next
            // wavelet.
            nWavelets++;
        }

        if((probe & 0x8000)==0x8000)
        {

            // anSamp is the number of waveform samples in the analog trace 
            unsigned short anSamp1 = buffer[0];
            evtfile.read((char*)buffer,BufferBytes);
            unsigned short anSamp2 = buffer[0];
            evtfile.read((char*)buffer,BufferBytes);
            const unsigned int anSamp = (anSamp2 << 16) | anSamp1;

            if((probe & 0x0003)==0x0002)
            {
                // analog probe is CFD
                if(anSamp > 0)
                {

                    for(int i=0;i<anSamp;i++)
                    {

                        evtfile.read((char*)buffer,BufferBytes);

                    }

                    // done with this event; increment the wavelet counter to get ready for the next
                    // wavelet.
                    nCFDs++;
                }
            }

            if((probe & 0x0003)==0x0001)
            {
                // analog probe is baseline
                if(anSamp > 0)
                {
                    
                    for(int i=0;i<anSamp;i++)
                    {

                        evtfile.read((char*)buffer,BufferBytes);

                    }

                    // done with this event; increment the wavelet counter to get ready for the next
                    // wavelet.
                    nBaselines++;
                }
            }
        }


        // Populate histograms with time data

        outCT->Fill(timetag);
        outFT->Fill(finetime);

        // Check for coincidence with another event (on a different channel) and fill the ROOT tree
        // if all coincident channels have been detected.

        if(channel == 6)
            {
                te.finetimeL = (finetime*2000)/1024; // 2000 ps per 1024 units;
                te.lgQL = lgQ;
            }

        else if(channel == 7)
            {
                te.finetimeR = (finetime*2000)/1024; // 2000 ps per 1024 units;
                te.lgQR = lgQ;
            }

        if(timetag == timeCoin /*|| timetag == timeCoin+1 || timetag == timeCoin-1*/)
        {
            //cout << timeCoin << " " << channel << endl;
            channelCoin++;

            if (channelCoin == COINCIDENCE_REQ)
            {
                te.timetag = timetag*2; // 2 ns per unit
                //cout << te.timetag << " " << te.finetimeL << " " << te.finetimeR << endl;

                if ((te.lgQL + te.lgQR) > 1800 && (te.lgQL + te.lgQR) < 2800)
                {
                    tree->Fill();
                }

                channelCoin = 1; // reset the coincindence counter for next event
            }

        }

        else
        {
            //cout << timeCoin << " " << channel << endl;
            timeCoin = timetag;
        }

    }

    else if(evtype==2)
    {

        unsigned short nSamp1 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        unsigned short nSamp2 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        unsigned long nSamp = (nSamp2 << 16) | nSamp1;

        for(int i=0;i<nSamp;i++)
        {

            evtfile.read((char*)buffer,BufferBytes);

        }

        // done with this event; increment the wavelet counter to get ready for the next
        // wavelet.
        nWaveforms++;
    }

    else
    {
        cout << "ERROR: unknown value for event type" << endl;
    }

}

void processRun(string evtname)
{
    outCT = new TH1S("outCT","outCT",1000000,0,1000000); // a ROOT plot to show coarse time (CT) of events
    outFT = new TH1S("outFT","outFT",1023,0,1023); // a ROOT plot to show fine time (FT) of events
    outFTLR = new TH2S("outFTLR","outFTLR",1023,0,1023,1023,0,1023); // a ROOT plot to show fine time (FT) of events
    outFTLR->SetMarkerStyle(20);

    // create text output to examine time differences from one event to the next
    ofstream timeDiff ("sorted/timeDiff.txt");

    // attempt to process the event file
    ifstream evtfile;

    evtfile.open(evtname,ios::binary);

    if (!evtfile)
    {
        cout << "Failed to open " << evtname << ". Please check that the file exists and is listed properly in runsToSort.txt" << endl;
        abort();
    }

    else
    {
        cout << evtname << " opened successfully. Start reading events..." << endl;

        evtfile.read((char*)buffer,BufferBytes);

        point = buffer;

        // start looping through the evtfile for events
        while(!evtfile.eof())
        {
            readHeader(evtfile);

            outCT->Fill(2*timetag);

            // print the time difference between adjacent events to timeDiff.txt
            timeDiff << (timetag - timetagP)*2 << endl;
            timetagP = timetag;

            unpack(evtfile); 

            // Event finished
            nE++;
        }

        // Input file finished
        cout << "Finished processing event file" << endl;
        cout << "Total events: " << nE << endl;

    }

    evtfile.close();
}

int main(int argc, char* argv[])
{

    // create a ROOT file
    TFile *file; 
    file = new TFile("finetimeCheck.root","RECREATE");
    file->cd();

    tree = new TTree("tree","");
    tree->SetAutoSave(0);
    tree->Branch("timetag",&te.timetag,"timetag/I:finetimeL/F:finetimeR/F");

    stringstream evtname;

    // open the event file(s)

    if (argc > 1) // list of runs given in the command line; IGNORE runsToSort.txt
    {
        for (int i=1; i<argc; i++)
        {
            evtname.str("");
            evtname << "../output/run" << argv[i] << ".evt";
            processRun(evtname.str());
        }
    }

    else
    {
        // open the event files
        ifstream evtFilenames;
        string evtFilename = "runsToSort.txt";
        ifstream evtfile;
        stringstream evtname;
        evtname << "../output/wutest.evt";
        string runNo = "-1";

        evtFilenames.clear();
        evtFilenames.open(evtFilename.c_str());

        if (!evtFilenames)
        {
            cout << "List of run numbers for sorting failed to open. Please check the input file (expected runsToSort.txt)" << endl;
            abort();
        }

        else // run number list is good - start processing runs.
        {
            while (getline(evtFilenames,runNo))
            {
                evtname.str("");
                evtname << "../output/run" <<  runNo << ".evt";
                processRun(evtname.str());
            }

        }
    }

    tree->Draw("finetimeL:finetimeR>>finetime2D","","SetMarkerStyle(20)");
    tree->Draw("finetimeR-finetimeL>>finetime1D","","SetMarkerStyle(20)");
    //tree->Scan("timetag:finetimeL:finetimeR");

    file->Write();
    tree->Write();
}
