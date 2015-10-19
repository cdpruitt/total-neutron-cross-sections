// Reads output from washudaq-X.X and sorts into ROOT and text files
// 
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
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TH1S.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectoryFile.h"
#include "MyClass.h"
#include <dirent.h>
#include <algorithm>
#include "sys/stat.h"
#include "unistd.h"
#include "time.h"
#include <regex>
#include <limits>

using namespace std;

// Buffer variables
int const BufferWords = 1;
int const BufferBytes = BufferWords*2;
unsigned short buffer[BufferWords];
unsigned short *point;

// Header fields for events
unsigned int size;
unsigned int evtType;
unsigned int channel;
unsigned int timetag;

// timetagP keeps track of the previous event's timetag, so we can count 
// macropulses by looking at timetag resets. Because channels don't read out in
// order, we need to track the most recent timetag for each enabled channel
vector<int> timetagP (8,0); // 8 channels all start with previous timetag = 0
vector<short> macroNo (8,0);  // 8 channels all start on the 0th macro
vector<short> evtNo (8,0);    // 8 channels all start on the 0th event

unsigned int nE = 0; // counter for the total number of events
unsigned int nWavelets = 0; // counter for the number of waveforms in DPP mode
unsigned int nCFDs = 0; // counter for the number of CFD traces (analog probe mode) 
unsigned int nBaselines = 0; // counter for the number of baseline traces (analog probe mode)
unsigned int nWaveforms = 0; // counter for the number of WAVEFORM mode waveforms

unsigned short sgQ, lgQ, fineTime, nSamp, probe, anSamp, extraSelect, extras1, extras2, puRej;
unsigned short chNo;

std::string runNo;

std::vector<short> waveform; // for holding one event's waveform data
std::vector<short> anProbe; // for holding one event's analog probe waveform data

std::vector<TH1S*> listWaveforms;
std::vector<TH1S*> listWavelets;
std::vector<TH1S*> listCFDs;
std::vector<TH1S*> listBaselines;

// for text file displaying events, organized by input channel to make sense of the data
ofstream totalOut ("sorted/allChannels.txt");
ofstream targetChangerOut ("sorted/targetChanger.txt");
ofstream monitorOut ("sorted/monitor.txt");
ofstream detectorLOut ("sorted/detectorL.txt");
ofstream detectorROut ("sorted/detectorR.txt");
ofstream detectorTOut ("sorted/detectorT.txt");

// for text output to examine time differences from one event to the next
ofstream timeDiff ("sorted/timeDiff.txt");

// A ROOT tree holds DPP data, keyed by a timestamp (timetag).
// Entries from different channels with the same timetag are associated via the tree.
TTree* tree;

struct event
{
    // label each event by runNo, macroNo, evtNo to uniquely identify
    unsigned short runNo, macroNo, evtNo;

    unsigned short chNo, evtType; // describe the event data

    unsigned int timetag;
    unsigned short fineTime, sgQ, lgQ;

    std::vector<short> waveform; // include the waveforms for each channel of the event
    std::vector<short> anProbe; // include the analog probes for each channel of the event
} ev;

bool text = false; // flag for producing text-file output apart from
                   // the default ROOT plots and ROOT tree
bool runlist = false; // flag indicating that runs should be read from
                      // runsToSort.txt, and NOT just the most recent file

TDirectoryFile *targetChangerDir, *monitorDir, *detectorLDir, *detectorRDir, *detectorTDir;

// Histograms for DPP data
TH1S* outTime; // coarse time
TH1S* outSGQ; // short gate Q
TH1S* outLGQ; // long gate Q
TH1S* outPZC; // positive zero-crossing for the CFD
TH1S* outNZC; // negative zero-crossing for the CFD
TH1S* outBaseline;
TH1S* outFT; // fine time
TH1S* outCT; // coarse time

// read the EVENT HEADER
// also determines where this event's data belongs (using channel)
int readHeader(ifstream& evtfile)
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

// prints EVENT HEADER to the appropriate text file
// 'out' points to a channel-specific text file
void printHeader(ofstream& out)
{
    stringstream temp; // used for formatting strings with units

    out << setfill('*') << setw(63) << "*" << endl;
    out << "| EVENT " << left << setfill(' ') << setw(54) << nE << "|" << endl;
    out << "|" << right << setfill('-') << setw(62) << "|" << endl;
    out << "| channel #" << channel;

    if(evtType==1)
    {
        out << left << setfill(' ') << setw(50) << ", DPP mode" << "|" << endl;
    }

    else if (evtType==2)
    {
        out << left << setfill(' ') << setw(50) << ", waveform mode " << "|" << endl;
    }

    else
    {
        out << left << setfill(' ') << setw(50) << ", ERROR DETERMINING MODE" << " |" << endl;
        cout << "Error: event type value out-of-range (DPP=1, waveform=2)" << endl;
        cout << "Event number = " << evtNo[chNo] << endl;
    }

    temp << size << " bytes";
    out << "| size = " << left << setfill(' ') << setw(53) << temp.str() << "|" << endl;

    temp.str("");
    temp << 2*timetag << " ns";
    out << "| timetag = " << left << setw(50) << temp.str() << "|" << endl;

    out << "|" << right << setfill('-') << setw(62) << "|" << endl;

    outCT->Fill(2*timetag);

}

// read the EVENT BODY
void readBody(ifstream& evtfile)
{
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

        // retrieve fine time from bits 0:9 (0x03ff)
        fineTime = (extras1 & 0x03ff);

        // sgQ is the short gate integrated charge, in digitizer units 
        sgQ = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);

        // lgQ is the long gate integrated charge, in digitizer units 
        lgQ = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        
        // baseline is the baseline level, frozen at the trigger time
/*        unsigned short baseline = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
*/

        // puRej is a pile-up rejection flag (1 if pile-up detected? 0 if not?)
        // Not reliable in current washudaq version - but still need to read the
        // word from the event file to avoid getting mismatched.
        puRej = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);

        // probe indicates whether an additional analog waveform will be captured along with the
        // input trace. The top bit turns on the analog probe; the bottom two bits describe the
        // probe type.
        probe = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        
        // nSamp is the number of waveform samples that follow (in LIST mode, this is 0)
        unsigned short nSamp1 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        unsigned short nSamp2 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        nSamp = (nSamp2 << 16) | nSamp1;

        waveform.clear();
        if(nSamp > 0)
        {
            // read waveform, consisting of nSamp samples
            for(int i=0;i<nSamp;i++)
            {
                waveform.push_back(buffer[0]);
                evtfile.read((char*)buffer,BufferBytes);
            }
        }

        if((probe & 0x8000)==0x8000) 
        {
            // analog probe enabled; read analog probe waveform
            unsigned short anSamp1 = buffer[0];
            evtfile.read((char*)buffer,BufferBytes);
            unsigned short anSamp2 = buffer[0];
            evtfile.read((char*)buffer,BufferBytes);
            anSamp = (anSamp2 << 16) | anSamp1;

            if(anSamp > 0)
            {
                for(int i=0;i<anSamp;i++)
                {
                    anProbe.push_back(buffer[0]);
                    evtfile.read((char*)buffer,BufferBytes);
                }
            }
        }
    }

    else if(evtType==2)
    {
        unsigned short nSamp1 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        unsigned short nSamp2 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        nSamp = (nSamp2 << 16) | nSamp1;

        waveform.clear();
        if(nSamp >0)
        {
            // read waveform, consisting of nSamp samples
            for(int i=0;i<nSamp;i++)
            {
                waveform.push_back(buffer[0]);
                evtfile.read((char*)buffer,BufferBytes);
            }
        }
    }

    else
    {
        cout << "ERROR: unknown value for event type" << endl;
        return;
    }
}

void printBody(ofstream& out)
{
    ostringstream histName;

    if(evtType==1)
    {
        stringstream temp;

        switch(extraSelect)
        {
            case 0:
                temp << extras2 << " *8.59 s";
                out << "| extended time stamp = " << left << setfill(' ') << setw(38) << temp.str() << "|" << endl;

                // readout gives baseline*4, and we want just baseline
                extras1 /= 4;
                out << "| baseline = " << left << setfill(' ') << setw(49) << extras1 << "|" << endl;
                break;

            case 1:
                temp << extras2 << " *8.59 s";
                out << "| extended time stamp = " << left << setfill(' ') << setw(34) << temp.str() << "|" << endl;
                // extract flags from bits 10:15 (0xfc00)
                out << "| flags = " << left << setfill(' ') << setw(49) << (extras1 & 0xfc00) << "|" << endl;
                break; 

            case 2:
                temp << extras2 << " *8.59 s";
                out << "| extended time stamp = " << left << setfill(' ') << setw(38) << temp.str() << "|" << endl;
                // extract flags from bits 10:15 (0xfc00)
                out << "| flags = " << left << setfill(' ') << setw(52) << (extras1 & 0xfc00) << "|" << endl;

                out << "| fine time stamp = " << left << setfill(' ') << setw(42) << fineTime << "|" << endl;
                outFT->Fill((extras1 & 0x03ff));
                break;

            case 3:
                // pulse peak value documentation
                out << "| pulse peak value" << left << setfill(' ') << setw(49) << extras1 << "|" << endl;
                break;

            case 5:
                // PZC and NZC
                out << "| PZC = " << left << setfill(' ') << setw(54) << extras2 << "|" << endl;
                out << "| NZC = " << left << setfill(' ') << setw(54) << extras1 << "|" << endl;
                outPZC->Fill(extras2);
                outNZC->Fill(extras1);
                break;

            case 7:
                // fixed value of 0x12345678
                const unsigned int extras = (extras2 << 16) | extras1;
                out << "| 305419896 ?= " << left << setfill(' ') << setw(47) << extras << "|" << endl;
                break;
        }

        out << "| short gate charge = " << left << setw(40) << sgQ << "|" << endl;
        out << "| long gate charge = " << left << setw(41) << lgQ << "|" << endl;

        // Pile-up detection not operational in current washudaq - ignore it.
        //out << "| pile-up detected = " << left << setw(40) << puRej << " |" << endl;

        temp.str("");
        temp << nSamp << " samples";
        out << "| waveform length = " << left << setw(41) << temp.str() << " |" << endl;
        out << "|" << right << setfill('-') << setw(62) << "|" << endl;

        if(nSamp > 0)
            // print wavelet data
        {
            out << left << setfill(' ') << setw(62) << "| Waveform samples" << "|" << endl;
            out << "|" << right << setfill(' ') << setw(62) << "|" << endl;

            for(std::vector<int>::size_type i = 0; i != waveform.size(); i++)
            {
                if(i%100 == 0 && i>0)
                {
                    out << "|" << right << setfill(' ') << setw(62) << "|" << endl;
                }

                if(i%10 == 0)
                {
                    temp.str("");
                    temp << "|";
                }

                temp << right << setfill(' ') << setw(6) << waveform[i];

                if(i%10==9 || i==nSamp-1)
                {
                    out << left << setw(62) << temp.str() << "|" << endl;
                } 

            }
            out << "|" << right << setfill('-') << setw(62) << "|" << endl;
        }

        if((probe & 0x8000)==0x8000)
        {
            // print analog probe data
            out << "| Analog probe enabled" << right << setfill(' ') << setw(41) << "|" << endl;

            temp.str("");
            temp << nSamp << " samples";
            out << "| waveform length = " << left << setw(41) << temp.str() << " |" << endl;
            out << "|" << right << setfill('-') << setw(62) << "|" << endl;

            if((probe & 0x0003)==0x0002)
            {
                // analog probe is CFD
                if(anSamp > 0)
                {
                    out << left << setfill(' ') << setw(62) << "| CFD samples" << "|" << endl;
                    out << "|" << right << setfill(' ') << setw(62) << "|" << endl;

                    for(std::vector<int>::size_type i = 0; i != anProbe.size(); i++)
                    {
                        if(i%100 == 0 && i>0)
                        {
                            out << "|" << right << setfill(' ') << setw(62) << "|" << endl;
                        }

                        if(i%10 == 0)
                        {
                            temp.str("");
                            temp << "|";
                        }

                        temp << right << setfill(' ') << setw(6) << anProbe[i];

                        if(i%10==9 || i==anSamp-1)
                        {
                            out << left << setw(62) << temp.str() << "|" << endl;
                        } 

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
                    histName.str("");
                    histName << "outBaseline" << evtNo[chNo];
                    string tempHist = histName.str();

                    listBaselines.push_back(new TH1S(tempHist.c_str(),tempHist.c_str(),nSamp,0,nSamp*2));

                    out << left << setfill(' ') << setw(62) << "| Baseline samples" << "|" << endl;
                    out << "|" << right << setfill(' ') << setw(62) << "|" << endl;

                    for(int i=0;i<anSamp;i++)
                    {
                        listBaselines[nBaselines]->SetBinContent(i,buffer[0]);

                        if(i%100 == 0 && i>0)
                        {
                            out << "|" << right << setfill(' ') << setw(62) << "|" << endl;
                        }

                        if(i%10 == 0)
                        {
                            temp.str("");
                            temp << "|";
                        }

                        temp << right << setfill(' ') << setw(6) << buffer[0];

                        if(i%10==9 || i==anSamp-1)
                        {
                            out << left << setw(62) << temp.str() << "|" << endl;
                        } 
                    }
                }
            }
        }

        else
        {
            out << "| Analog probe disabled" << right << setfill(' ') << setw(40) << "|" << endl;
            out << "|" << right << setfill(' ') << setw(62) << "|" << endl;
            out << setfill('*') << setw(63) << "*" << endl;
        }
    }

    else if(evtType==2)
    {
        stringstream temp;
        temp << nSamp << " samples";
        out << "| waveform length = " << left << setfill(' ') << setw(41) << temp.str() << " |" << endl;
        out << "|" << right << setfill('-') << setw(62) << "|" << endl;

        for(int i=0;i<nSamp;i++)
        {

            if(nSamp>1000 && (i%1000 == 0))
            {
                out << "|" << right << setfill(' ') << setw(62) << "|" << endl;
                stringstream temp;
                temp << "Samples " << i << "-" << i+1000;
                out << "| " << left << setw(60) << temp.str() << "|" << endl;
            }

            if(i%100 == 0)
            {
                out << "|" << right << setfill(' ') << setw(62) << "|" << endl;
            }

            if(i%10 == 0)
            {
                out << "|";
            } 

            out << right << setfill(' ') << setw(6) << waveform[i];

            if(i%10 == 9)
            {
                out << " |" << endl;
            } 
        }

        out << "|" << right << setfill(' ') << setw(62) << "|" << endl;
    }

    else
    {
        cout << "ERROR: unknown value for event type" << endl;
        return;
    }
    out << endl;
}

void fillTree()
{
    if(timetag <= timetagP[chNo])
    {
        // new macropulse for this channel
        macroNo[chNo]++;
    }
    timetagP[chNo] = timetag;

    // detectors are summed; simply fill tree with the event
    ev.runNo = std::stoi(runNo);
    ev.macroNo = macroNo[chNo];
    ev.evtNo = evtNo[chNo];
    ev.chNo = chNo;
    ev.evtType = evtType;
    ev.timetag = timetag;
    ev.fineTime = fineTime;
    ev.sgQ = sgQ;
    ev.lgQ = lgQ;
    //ev.waveform = waveform;
    //ev.anProbe = anProbe;

    tree->Fill();
}

void processRun(string evtname)
{
    // attempt to process the event file
    ifstream evtfile;
    evtfile.open(evtname,ios::binary);

    if (!evtfile)
    {
        cout << "Failed to open " << evtname << ". Please check that the file exists" << endl;
        abort();
    }

    else // individual run is good - start processing events from that run.
    {
        cout << evtname << " opened successfully. Start reading events..." << endl;

        evtfile.read((char*)buffer,BufferBytes);

        point = buffer;
        std::fill(timetagP.begin(), timetagP.end(), 0); // reset timetagP so macroNo counting works properly
        std::fill(macroNo.begin(), macroNo.end(), 0); // reset macroNo, so each macro in a run is counted right

        // start looping through the evtfile for events
        while(!evtfile.eof())
        {
            // get channel number from the event header
            chNo = readHeader(evtfile);

            readBody(evtfile); // extract data from the event body
            fillTree(); // fill the tree with data from the event body

            if(text)
            {
                // text output enabled

                // print the header to an allChannels.txt
                printHeader(totalOut);

                // print the time difference between adjacent events to timeDiff.txt
                timeDiff << (timetag - timetagP[chNo])*2 << endl;
                timetagP[chNo] = timetag;

                switch (chNo)
                {
                    case 0: 
                        // target-changer data
                        targetChangerDir->cd();
                        printHeader(targetChangerOut);
                        printBody(targetChangerOut);
                        break;

                    case 1:
                        break;

                    case 2:
                        // monitor data
                        monitorDir->cd();
                        printHeader(monitorOut);
                        printBody(monitorOut);
                        break;

                    case 3:
                        break;

                    case 4:
                        // Detector data (assume detectors T'd together)
                        detectorTDir->cd();
                        printHeader(detectorTOut);
                        printBody(detectorTOut);
                        break;

                    case 5:
                        break;

                    case 6:
                        // Detector data (left detector only)
                        detectorLDir->cd();
                        printHeader(detectorLOut);
                        printBody(detectorLOut);
                        break;

                    case 7:
                        // Detector data (right detector only)
                        detectorRDir->cd();
                        printHeader(detectorROut);
                        printBody(detectorROut);
                        break;

                    default:
                        cout << "ERROR: unknown value for channel type" << endl;
                        break;
                }
            }

            // Event finished
            evtNo[chNo]++;
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

    // create a ROOT file for holding events, with one directory per input channel
    TFile *file; 
    file = new TFile("events.root","RECREATE");
    file->cd();

    tree = new TTree("tree","");
    tree->SetAutoSave(0);
    tree->Branch("runNo",&ev.runNo,"runNo/s");
    tree->Branch("macroNo",&ev.macroNo,"macroNo/s");
    tree->Branch("evtNo",&ev.evtNo,"evtNo/s");
    tree->Branch("chNo",&ev.chNo,"chNo/s");
    tree->Branch("evtType",&ev.evtType,"evtType/i");
    tree->Branch("timetag",&ev.timetag,"timetag/i");
    tree->Branch("fineTime",&ev.fineTime,"fineTime/s");
    tree->Branch("sgQ",&ev.sgQ,"sgQ/s");
    tree->Branch("lgQ",&ev.lgQ,"lgQ/s");

    stringstream evtname;
    
    // define ROOT histograms
    outTime = new TH1S("outTime","outTime",10000,0,10000000000);
    outSGQ = new TH1S("outSGQ","outSGQ",1024,0,70000);
    outLGQ = new TH1S("outLGQ","outLGQ",1024,0,70000);
    outPZC = new TH1S("outPZC","outPZC",16384,0,16384);
    outNZC = new TH1S("outNZC","outNZC",16384,0,16384);
    outBaseline = new TH1S("outBaseline","outBaseline",1024,0,17000);
    outFT = new TH1S("outFT","outFT",1023,0,1023); // a ROOT plot to show fine time (FT) of events
    outCT = new TH1S("outCT","outCT",1000000,0,1000000); // a ROOT plot to show coarse time (CT) of events

    // set ROOT directory structure
    targetChangerDir = new TDirectoryFile("targetChanger","Target Changer");
    monitorDir = new TDirectoryFile("monitor","Monitor");
    detectorLDir = new TDirectoryFile("detL","Detector left)");
    detectorRDir = new TDirectoryFile("detR","Detector right");
    detectorTDir = new TDirectoryFile("detT","Detector sum");
    // TDirectory *dirSub; to be uncommented and used for sub-directories, if desired

    // open the event files

    if (argc > 1) // flags detected; modify behavior based on flags
    {
        for (int i=1; i<argc; i++)
        {
            if (std::string(argv[i]) == "--text" || std::string(argv[i]) == "-t")
            {
                // produce text files for each channel containing all event
                // data from the input file. This will significantly increase
                // processing time for this program
                text = true;
            }

            if (std::string(argv[i]) == "--runlist" || std::string(argv[i]) == "-rl")
            {
                // read the runs to be processed from runsToSort.txt instead
                // of using the most recently modified file in ../output
                runlist = true;
            }

        }
    }

    if (runlist)
    {
        ifstream evtFilenames;
        string evtFilename = "runsToSort.txt";

        evtFilenames.clear();
        evtFilenames.open(evtFilename.c_str());

        if (!evtFilenames)
        {
            cout << "The list of event files failed to open. Please check the input file (runsToSort.txt)" << endl;
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

    else
    {
        FILE *fp;
        char path[100];

        // Open the command for reading files
        fp = popen("ls -t ../output | head -1 | egrep -o '[0-9]+'", "r");
        if (fp == NULL)
        {
            std::cout  << "Failed to find most recent file in ../output" << std::endl;
            exit(1);
        }

        fscanf(fp,"%s",path);
        evtname << "../output/run" << path << ".evt";

        std::stringstream temp;
        temp << path;
        runNo = temp.str();

        processRun(evtname.str());
    }

    file->Write();
}
