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
#include "TFile.h"
#include "TTree.h"
#include "TEntryList.h"
#include "TH1I.h"
#include "TDirectoryFile.h"
#include <dirent.h>
#include <algorithm>
#include "sys/stat.h"
#include "unistd.h"
#include "time.h"
#include <regex>
#include <limits>
#include "TROOT.h"
#include "TRandom3.h"

using namespace std;

// Experimental constants


const float FLIGHT_DISTANCE = 2500; // detector distance from source, in cm

// Target positions and identities
const string TARGETS[6] = {"blank","carbonS","carbonL","112Sn","Nat Sn","124Sn"};

// Time delay of target changer coarse time after real macropulse start time
const double TIME_OFFSET = 940; // in ns

// Period of micropulses
const double MICRO_PERIOD = 1788.82; // in ns

// Target Changer lgQ gates for setting integral target positions 1-6
const int tarGate[12] = {5000,7500,11000,13500,17000,20000,21000,25000,27000,31000,32000,36000};

// Number of wavelets in waveform mode that constitute a complete macropulse
const int WAVELET_NO = 11;

// Physical constants
const float C = 299792458; // speed of light in m/s

const float NEUTRON_MASS = 939.56536; // in MeV/c^2



// Buffer variables
int const BufferWords = 1;
int const BufferBytes = BufferWords*2;
unsigned short buffer[BufferWords];
unsigned short *point;

// Header fields for events
unsigned int size;
unsigned int evtType;
unsigned int chNo;
unsigned int timetag;
unsigned int extTime;

// timetagP keeps track of the previous event's timetag, so we can count 
// macropulses by looking at timetag resets. Because channels don't read out in
// order, we need to track the most recent timetag for each enabled channel
vector<int> timetagP (8,0); // 8 channels all start with previous timetag = 0
vector<int> evtNo (8,0);    // 8 channels all start on the 0th event

// used to perform special behavior when the first wavelet of a new waveform
// is detected during reconstruction of a full macropulse waveform

unsigned int nE = 0; // counter for the total number of events
unsigned int macroNo = 0; // keep track of target changer macropulse number
unsigned int tempMacro = 0; // keep track of an event's macropulse number

/*unsigned int nWavelets = 0; // counter for the number of waveforms in DPP mode
unsigned int nCFDs = 0; // counter for the number of CFD traces (analog probe mode) 
unsigned int nBaselines = 0; // counter for the number of baseline traces (analog probe mode)
unsigned int nWaveforms = 0; // counter for the number of WAVEFORM mode waveforms
*/

std::string runNo;
std::string runDir;

unsigned int sgQ, lgQ, fineTime, nSamp, probe, anSamp, extraSelect, extras1, extras2, puRej;

std::vector<int> waveform; // for holding one event's waveform data
std::vector<int> anProbe; // for holding one event's analog probe waveform data

// for text file displaying events, organized by input channel to make sense of the data
ofstream totalOut ("textSort/allChannels.txt");
ofstream targetChangerOut ("textSort/targetChanger.txt");
ofstream monitorOut ("textSort/monitor.txt");
ofstream detectorLOut ("textSort/detectorL.txt");
ofstream detectorROut ("textSort/detectorR.txt");
ofstream detectorTOut ("textSort/detectorT.txt");

vector<TDirectoryFile*> directs;
TH1I* DPPWaveform;
TH1I* WaveWaveform0;
TH1I* WaveWaveform2;
TH1I* WaveWaveform4;
TH1I* WaveWaveform6;
TH1I* WaveWaveform7;
int waveformStart[8]; // for indicating the timetag of the first waveform

// ROOT file directory structure 
string dirs[8] = {"targetChanger","","monitor","","detT","","detL","detR"};
TDirectory *DPPWaveformsDir;
TDirectory *WaveWaveformsDir;

vector<vector<TH1I*>> histos; // holds all histograms for the run
// split into sub-vectors on a per-channel basis

TH1I* outMacro;
TH1I* outEvt;
TH1I* outExtTime;
TH1I* outTime;
TH1I* outFT;
TH1I* outSGQ;
TH1I* outLGQ;

// for text output to examine time differences from one event to the next
ofstream timeDiff ("textSort/timeDiff.txt");

// Data is sorted into a ROOT tree; each event has a unique runNo-macroNo-evtNo ID.
TTree* tree;

TRandom3* rng;

struct event
{
    // label each event by runNo, macroNo, evtNo to uniquely identify
    unsigned int runNo, macroNo, evtNo;

    unsigned int chNo; // describe data stream origin (i.e., detector)
    unsigned int evtType; // describe the event data: either DPP or waveform

    unsigned int timetag;
    unsigned int extTime;
    unsigned int fineTime, sgQ, lgQ;

    vector<int> waveform; // include the waveforms for each channel of the event
    vector<int> anProbe; // include the analog probes for each channel of the event
} ev;

// hold previous 500 target times and macro numbers
vector<pair<unsigned int,double>> targetTimeList;
double targetTime; // this is the correct target time for each detector event
double prevTime = 0;
int prevTarget = 1;

vector<bool> firstWaveform (8,false);

//TEventList *targetCh; // for sorting events by their target positions

bool text = false; // flag for producing text-file output apart from
                   // the default ROOT plots and tree filling
bool cs = false; // flag indicating that a cross-section plot should be
                 // produced for each target position

TDirectoryFile *targetChangerDir, *monitorDir, *detectorLDir, *detectorRDir, *detectorTDir;

// read the EVENT HEADER
// also determines where this event's data belongs (using chNo)
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

    //timetagP[chNo] = timetag;

    return chNo;
}

// prints EVENT HEADER to the appropriate text file
// 'out' points to a channel-specific text file
void printHeader(ofstream& out)
{
    stringstream temp; // used for formatting strings with units

    out << setfill('*') << setw(63) << "*" << endl;
    out << "| EVENT " << left << setfill(' ') << setw(54) << nE << "|" << endl;
    out << "|" << right << setfill('-') << setw(62) << "|" << endl;
    out << "| run " << runNo << ", macro " << left << setfill(' ') << setw(44) << macroNo << "|" << endl;
    out << "| channel " << chNo;

    if(evtType==1)
    {
        out << left << setfill(' ') << setw(51) << ", DPP mode" << "|" << endl;
    }

    else if (evtType==2)
    {
        out << left << setfill(' ') << setw(51) << ", waveform mode " << "|" << endl;
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
    temp << timetag << " ns";
    out << "| timetag = " << left << setw(50) << temp.str() << "|" << endl;

    out << "|" << right << setfill('-') << setw(62) << "|" << endl;


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

        switch(extraSelect)
        {
            case 2:
                // retrieve extended time from bits 16-31 (extras2)
                extTime = extras2;
                // retrieve fine time from bits 0:9 (0x03ff)
                fineTime = (extras1 & 0x03ff);
                break;

            default:
                ;
                // unimplemented
        }

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
                break;

            case 3:
                // pulse peak value documentation
                out << "| pulse peak value" << left << setfill(' ') << setw(49) << extras1 << "|" << endl;
                break;

            case 5:
                // PZC and NZC
                out << "| PZC = " << left << setfill(' ') << setw(54) << extras2 << "|" << endl;
                out << "| NZC = " << left << setfill(' ') << setw(54) << extras1 << "|" << endl;
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
                    //nCFDs++;
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

                    //listBaselines.push_back(new TH1I(tempHist.c_str(),tempHist.c_str(),nSamp,0,nSamp*2));

                    out << left << setfill(' ') << setw(62) << "| Baseline samples" << "|" << endl;
                    out << "|" << right << setfill(' ') << setw(62) << "|" << endl;

                    for(int i=0;i<anSamp;i++)
                    {
                        //listBaselines[nBaselines]->SetBinContent(i,buffer[0]);

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
    //ev.runNo = std::stoi(runNo);
    //ev.macroNo = macroNo;
    //ev.evtNo = evtNo[chNo];
    ev.chNo = chNo;
    ev.evtType = evtType;
    ev.extTime = extTime;
    ev.timetag = timetag;
    ev.fineTime = fineTime;
    ev.sgQ = sgQ;
    ev.lgQ = lgQ;
    //ev.waveform = waveform;
    //ev.anProbe = anProbe;

    tree->Fill();

}

void fillHistos()
{
    int totalEntries = tree->GetEntries();

    targetTime = targetTimeList.front().second;

    for (int i=0; i<totalEntries; i++)
    {
        tree->GetEntry(i);

        gDirectory->cd("/");
        gDirectory->GetDirectory(dirs[chNo].c_str())->cd();

        // don't plot events coming >1ms after a target changer event
        if ((timetag+pow(2,32)*extTime)-targetTime < 1000000) 
        {
            if ((timetag+pow(2,32)*extTime)-targetTime < 0) 
            {
                cout << timetag+pow(2,32)*extTime << " " << targetTime << endl;
            }

            if (evtType == 1)
            {
                // DPP mode - fill DPP histograms

                outMacro = (TH1I*)(gDirectory->FindObject("outMacro"));
                outMacro->Fill(tempMacro);

                outEvt = (TH1I*)(gDirectory->FindObject("outEvt"));
                outEvt->Fill(evtNo[chNo]);

                outExtTime = (TH1I*)(gDirectory->FindObject("outExtTime"));
                outExtTime->Fill(extTime);

                outTime = (TH1I*)(gDirectory->FindObject("outTime"));
                outTime->Fill(timetag);

                outSGQ = (TH1I*)(gDirectory->FindObject("outSGQ"));
                outSGQ->Fill(sgQ);

                outLGQ = (TH1I*)(gDirectory->FindObject("outLGQ"));
                outLGQ->Fill(lgQ);

                outFT = (TH1I*)(gDirectory->FindObject("outFT"));
                outFT->Fill(fineTime + rng->Rndm());

                stringstream temp;
                temp << "macroNo " << macroNo << "evtNo " << evtNo[chNo];

                DPPWaveformsDir = (TDirectory*)gDirectory->FindObject("DPPWaveformsDir");
                DPPWaveformsDir->cd();

                if (evtNo[chNo]%10000 == 0)
                {
                    DPPWaveform = new TH1I(temp.str().c_str(),temp.str().c_str(),waveform.size(),0,waveform.size()*2);
                    //cout << timetag+pow(2,32)*extTime-targetTime << endl;

                    for(int i=0;i<waveform.size();i++)
                    {
                        DPPWaveform->SetBinContent(i,waveform[i]);
                    }
                }

            }

            else if (evtType == 2)
            {
                // waveform mode - create full macropulse waveforms

                WaveWaveformsDir = (TDirectory*)gDirectory->FindObject("WaveWaveformsDir");
                WaveWaveformsDir->cd();

                stringstream temp;

                switch (chNo)
                {
                    case 0:

                        if (timetag > waveformStart[chNo] + 1000000)
                        {
                            // new waveform-mode macropulse
                            temp << "macropulse " << macroNo;

                            WaveWaveform0 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);

                            waveformStart[chNo] = timetag;
                        }

                        for(int i=0;i<waveform.size();i++)
                        {
                            WaveWaveform0->SetBinContent(i+(timetag-waveformStart[chNo])/2,waveform[i]);
                        }

                        break;

                    case 2:

                        if (timetag > waveformStart[chNo] + 1000000)
                        {
                            temp << "macropulse " << macroNo;

                            WaveWaveform2 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);

                            waveformStart[chNo] = timetag;
                        }


                        for(int i=0;i<waveform.size();i++)
                        {
                            WaveWaveform2->SetBinContent(i+(timetag-waveformStart[chNo])/2,waveform[i]);
                        }

                        break;

                    case 4:

                        if (timetag > waveformStart[chNo] + 1000000)
                        {
                            temp << "macropulse " << macroNo;

                            WaveWaveform4 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);

                            waveformStart[chNo] = timetag;
                        }


                        for(int i=0;i<waveform.size();i++)
                        {
                            WaveWaveform4->SetBinContent(i+(timetag-waveformStart[chNo])/2,waveform[i]);
                        }

                        break;

                    case 6:

                        if (timetag > waveformStart[chNo] + 1000000)
                        {

                            temp << "macropulse " << macroNo;

                            WaveWaveform6 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);

                            waveformStart[chNo] = timetag;
                        }


                        for(int i=0;i<waveform.size();i++)
                        {
                            WaveWaveform6->SetBinContent(i+(timetag-waveformStart[chNo])/2,waveform[i]);
                        }

                        break;

                    case 7:

                        if (timetag > waveformStart[chNo] + 1000000)
                        {
                            temp << "macropulse " << macroNo;

                            WaveWaveform7 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);

                            waveformStart[chNo] = timetag;
                        }


                        for(int i=0;i<waveform.size();i++)
                        {
                            WaveWaveform7->SetBinContent(i+(timetag-waveformStart[chNo])/2,waveform[i]);
                        }

                        break;
                }
            }
        }
}

void calculateCS()
{
    // create cross-section histograms

    fill(targetTimeList.begin(),targetTimeList.end(),make_pair(0,0)); // clear the target time list

    gDirectory->cd("/");

    TH1I *blankRaw = new TH1I("blank","blank",700000,0,700);
    TH1I *carbonSRaw = new TH1I("carbonS","carbonS",700000,0,700);
    TH1I *carbonLRaw = new TH1I("carbonL","carbonL",700000,0,700);
    TH1I *Sn112Raw = new TH1I("Sn112","Sn112",700000,0,700);
    TH1I *NatSnRaw = new TH1I("NatSn","NatSn",700000,0,700);
    TH1I *Sn124Raw = new TH1I("Sn124","Sn124",700000,0,700);
    TH1I *totalRaw = new TH1I("total","total",700000,0,700);

    TH1I *TTOF = new TH1I("TTOF","All events time of flight",100000,0,2000);
    TH1I *MTOF = new TH1I("MTOF","Monitor time of flight",100000,0,2000);
    TH1I *STOF = new TH1I("STOF","Summed-detector time of flight",100000,0,2000);
    TH1I *LTOF = new TH1I("LTOF","Left detector time of flight",100000,0,2000);
    TH1I *RTOF = new TH1I("RTOF","Right detector time of flight",100000,0,2000);

    // link tree branches to variables-to-read
    //tree->SetBranchAddress("macroNo",&macroNo);
    tree->SetBranchAddress("chNo",&chNo);
    tree->SetBranchAddress("lgQ",&lgQ);
    tree->SetBranchAddress("timetag",&timetag);
    tree->SetBranchAddress("fineTime",&fineTime);
    tree->SetBranchAddress("extTime",&extTime);

    int targetPos = -1;

    int totalEntries = tree->GetEntries();
    
    for (int i=0; i<totalEntries; i++)
    {
        tree->GetEntry(i);

        if (evtType==1) // don't use waveform events to populate
        {

            if (chNo == 0) // update target position and target changer time
            {
                // assign an integral target position based on lgQ of target changer signal
                if (lgQ>tarGate[0] && lgQ<tarGate[1])
                {
                    targetPos = 1;
                }

                else if (lgQ>tarGate[2] && lgQ<tarGate[3])
                {
                    targetPos = 2;
                }

                else if (lgQ>tarGate[4] && lgQ<tarGate[5])
                {
                    targetPos = 3;
                }

                else if (lgQ>tarGate[6] && lgQ<tarGate[7])
                {
                    targetPos = 4;
                }

                else if (lgQ>tarGate[8] && lgQ<tarGate[9])
                {
                    targetPos = 5;
                }

                else if (lgQ>tarGate[10] && lgQ<tarGate[11])
                {
                    targetPos = 6;
                }

                else
                {
                    targetPos = -1;
                }

            }

            else if (chNo == 2 || chNo == 4 || chNo == 6 || chNo == 7) // detector event; continue
            {
                // calculate a targetTime-calibrated time-of-flight of neutrons
                // trueTime is neutron time-of-flight since micropulse
                // start, in ns

                for(vector<pair<unsigned int,double>>::reverse_iterator it = targetTimeList.rbegin(); it != targetTimeList.rend(); ++it)
                {
                    if ((timetag+pow(2,32)*extTime)>it->second && (timetag+pow(2,32)*extTime)<it->second+1000000) 
                    {
                        // found a previous target time within 1 ms 
                        targetTime = it->second;
                        break;
                    }
                }

                float trueTime = fmod(((double)timetag+((double)fineTime*2./1024.)-targetTime-TIME_OFFSET),MICRO_PERIOD);

                // convert trueTime into neutron velocity based on flight path distance
                float velocity = pow(10,7)*FLIGHT_DISTANCE/trueTime; // in meters/sec 

                // convert velocity to relativistic kinetic energy
                float rKE = (pow((1-pow((velocity/C),2)),-0.5)-1)*NEUTRON_MASS; // in MeV

                if (trueTime > 100) // gate disallowing gammas
                {
                    switch (targetPos)
                    {
                        case 1:
                            // BLANK
                            blankRaw->Fill(rKE);
                            break;
                        case 2:
                            // SHORT CARBON
                            carbonSRaw->Fill(rKE);
                            break;
                        case 3:
                            // LONG CARBON
                            carbonLRaw->Fill(rKE);
                            break;
                        case 4:
                            // Sn112
                            Sn112Raw->Fill(rKE);
                            break;
                        case 5:
                            // Natural Sn
                            NatSnRaw->Fill(rKE);
                            break;
                        case 6:
                            // Sn124
                            Sn124Raw->Fill(rKE);
                            break;
                        default:
                            break;
                    }
                    totalRaw->Fill(rKE);
                }

                TTOF->Fill(trueTime);

                switch (chNo)
                {
                    case 2:
                        MTOF->Fill(trueTime);
                        break;
                    case 4:
                        STOF->Fill(trueTime);
                        break;
                    case 6:
                        LTOF->Fill(trueTime);
                        break;
                    case 8:
                        RTOF->Fill(trueTime);
                }

                if (i%10000 == 0)
                {
                    //cout << "trueTime " << trueTime << " velocity " << velocity << " rKE " << rKE << endl;
                    cout << "Populated " << i << " events in cross-section histograms\r";
                    fflush(stdout);
                }
            }
        }
    }
    cout << endl << endl;

    //TH1I* carbonS = new TH1I("
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
        //std::fill(timetagP.begin(), timetagP.end(), 0); // reset timetagP so macroNo counting works properly
        //std::fill(macroNo.begin(), macroNo.end(), 0); // reset macroNo, so each macro in a run is counted right

        // start looping through the evtfile for events

        int looplimit = 0;

        while(!evtfile.eof() && looplimit<10000000)
        {
            looplimit++;
            // get channel number from the event header
            chNo = readHeader(evtfile);

            readBody(evtfile); // extract data from the event body

            fillTree(); // fill the tree with event data

            if(text)
            {
                // text output enabled

                // print the header to an allChannels.txt
                printHeader(totalOut);

                // print the time difference between adjacent events to timeDiff.txt
                //timeDiff << (timetag - timetagP[chNo])*2 << endl;
                //timetagP[chNo] = timetag;

                switch (chNo)
                {
                    case 0: 
                        // target-changer data
                        //targetChangerDir->cd();
                        printHeader(targetChangerOut);
                        printBody(targetChangerOut);
                        break;

                    case 1:
                        break;

                    case 2:
                        // monitor data
                        //monitorDir->cd();
                        printHeader(monitorOut);
                        printBody(monitorOut);
                        break;

                    case 3:
                        break;

                    case 4:
                        // Detector data (assume detectors T'd together)
                        //detectorTDir->cd();
                        printHeader(detectorTOut);
                        printBody(detectorTOut);
                        break;

                    case 5:
                        break;

                    case 6:
                        // Detector data (left detector only)
                        //detectorLDir->cd();
                        printHeader(detectorLOut);
                        printBody(detectorLOut);
                        break;

                    case 7:
                        // Detector data (right detector only)
                        //detectorRDir->cd();
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

void populateMacros()
{
    int totalEntries = tree->GetEntries();

    for (int i=0; i<totalEntries; i++)
    {
        tree->GetEntry(i);

        // Refresh targetCh time
        if (chNo==0) 
        {
            if (evtType==1)
            {
                // new macropulse in DPP mode
                macroNo++;
                fill_n(evtNo.begin(),8,0); // clear the evtNo counter for each channel
            }

            else if (evtType==2)
            {
                if (prevTarget==1)
                {
                    macroNo++;
                    fill_n(evtNo.begin(),8,0); // clear the evtNo counter for each channel

                    prevTime = timetag+pow(2,32)*extTime;
                }

                else if (timetag+pow(2,32)*extTime>prevTime+7000000)
                {
                    macroNo++;
                    fill_n(evtNo.begin(),8,0); // clear the evtNo counter for each channel
                }

            }

            targetTimeList.push_back(make_pair(macroNo,timetag + pow(2,32)*extTime));
            prevTarget = evtType;
        }
    }
}

int main(int argc, char* argv[])
{

    /*************************************************************************/
    /* read flags to set the sorting mode                                    */
    /*************************************************************************/ 

    if (argc > 3) // flags detected
    {
        if (std::string(argv[3]) == "true")
        {
            // produce text files for each channel containing all event
            // data from the input file. This will significantly increase
            // processing time for this program
            text = true;
        }

        if (string(argv[4]) == "true")
        {
            // produce cross-section plots based on target position
            cs = true;
        }
    }

    runDir = argv[1];
    runNo = argv[2];

    // Create a tree for this run
    TFile *file;

    stringstream treeName;
    stringstream fileName;
    treeName << runDir << "-" << runNo; 
    fileName << "/media/ExternalDrive1/analysis/" << runDir << "/" << treeName.str() << ".root";

    file = new TFile(fileName.str().c_str(),"UPDATE");

    if(file->Get("tree"))
    {
        tree = (TTree*)file->Get("tree");
        cout << "Found previous tree; skipping populating events into tree" << endl;
    }

    else
    {
        tree = new TTree("tree","");
        cout << "Created ROOT tree " << treeName.str() << endl;

        //tree->Branch("runNo",&ev.runNo,"runNo/i");
        //tree->Branch("macroNo",&ev.macroNo,"macroNo/i");
        //tree->Branch("evtNo",&ev.evtNo,"evtNo/i");
        tree->Branch("chNo",&ev.chNo,"chNo/i");
        tree->Branch("evtType",&ev.evtType,"evtType/i");
        tree->Branch("extTime",&ev.extTime,"extTime/i");
        tree->Branch("timetag",&ev.timetag,"timetag/i");
        tree->Branch("fineTime",&ev.fineTime,"fineTime/i");
        tree->Branch("sgQ",&ev.sgQ,"sgQ/i");
        tree->Branch("lgQ",&ev.lgQ,"lgQ/i");
        //tree->Branch("waveform",&ev.waveform);

        //ENABLE for diagnostic analog probe
        //tree->Branch("anProbe",&ev.anProbe);

        for(int i=0; i<8; i++)
        {
            vector<TH1I*> tempVec;
            histos.push_back(tempVec); // create sub-vector for this channel

            if (dirs[i].compare("") != 0) // valid channel
            {
                gDirectory->mkdir(dirs[i].c_str(),dirs[i].c_str());
                gDirectory->GetDirectory(dirs[i].c_str())->cd();

                // instantiate histograms

                histos.back().push_back(new TH1I("outMacro","outMacro",100000,0,100000));
                histos.back().back()->GetXaxis()->SetTitle("macropulse number");

                histos.back().push_back(new TH1I("outEvt","outEvt",1000,0,1000));
                histos.back().push_back(new TH1I("outExtTime","outExtTime",1000,0,1000));
                histos.back().push_back(new TH1I("outTime","outTime",250000,0,2500000000));
                histos.back().push_back(new TH1I("outSGQ","outSGQ",35000,0,35000));
                histos.back().push_back(new TH1I("outLGQ","outLGQ",70000,0,70000));
                histos.back().push_back(new TH1I("outFT","outFT",1023,0,1023));


                gDirectory->mkdir("DPPWaveformsDir","raw DPP waveforms");
                gDirectory->mkdir("WaveWaveformsDir","concatenated waveform waveforms");

                gDirectory->cd("/");
            }
        }

        rng = new TRandom3();

        stringstream runName;
        runName << "/media/ExternalDrive1/output/" << runDir << "/data-" << runNo << ".evt";
        processRun(runName.str());

        populateMacros();
        fillHistos(); // fill the histos with event data

    }

    if (cs)
    {
        calculateCS(); // produce histograms showing cross-sections for each target
    }

    file->Write();

}
