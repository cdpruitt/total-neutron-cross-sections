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

using namespace std;

// Experimental constants


const float FLIGHT_DISTANCE = 2500; // detector distance from source, in cm

// Target positions and identities
const string TARGETS[6] = {"blank","carbonS","carbonL","112Sn","Nat Sn","124Sn"};

// Time delay of target changer coarse time after real macropulse start time
<<<<<<< Updated upstream
const float TIME_OFFSET = 300; // in ns
=======
const double TIME_OFFSET = 951; // in ns
>>>>>>> Stashed changes

// Period of micropulses
const float MICRO_PERIOD = 1788.82; // in ns

// Target Changer lgQ gates for setting integral target positions 1-6
const int tarGate[12] = {5000,7500,10500,13500,16000,19000,21000,25000,27000,31000,32000,36000};

// Number of wavelets in waveform mode that constitute a complete macropulse
const int WAVELET_NO = 11;

// Physical constants
const float C = 299792458; // speed of light in m/s

const float NEUTRON_MASS = 939.56536; // in MeV/c^2

// active channels
int chNoList[5] = {0,2,4,6,7};


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

vector<bool> firstWaveform (8,false);
// used to perform special behavior when the first wavelet of a new waveform
// is detected during reconstruction of a full macropulse waveform

unsigned int nE = 0; // counter for the total number of events
<<<<<<< Updated upstream
unsigned int macroNo = 0;
=======
unsigned int tempMacro = 0; // keep track of an event's macropulse number
>>>>>>> Stashed changes

/*unsigned int nWavelets = 0; // counter for the number of waveforms in DPP mode
unsigned int nCFDs = 0; // counter for the number of CFD traces (analog probe mode) 
unsigned int nBaselines = 0; // counter for the number of baseline traces (analog probe mode)
unsigned int nWaveforms = 0; // counter for the number of WAVEFORM mode waveforms
*/

std::string runNo;
unsigned int sgQ, lgQ, fineTime, nSamp, probe, anSamp, extraSelect, extras1, extras2, puRej;

vector<int> waveform; // for holding one event's waveform data
vector<int> anProbe; // for holding one event's analog probe waveform data

vector<TDirectoryFile*> directs;
TH1I* DPPWaveform;
TH1I* WaveWaveform0;
TH1I* WaveWaveform2;
TH1I* WaveWaveform4;
TH1I* WaveWaveform6;
TH1I* WaveWaveform7;
vector<int> waveformStart (8,0); // for indicating the timetag of the first waveform

// ROOT file directory structure 
string dirs[5] = {"targetChanger","monitor","detT","detL","detR"};
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

TTree* tree;

struct event
{
    unsigned int chNo; // describe data stream origin (i.e., detector)
    unsigned int evtType; // describe the event data: either DPP or waveform

    unsigned int timetag, extTime, fineTime; // event timestamps
    unsigned int sgQ, lgQ; // event charge gates

    vector<int> waveform; // include the waveforms for each channel of the event
    vector<int> anProbe; // include the analog probes for each channel of the event
} ev;

<<<<<<< Updated upstream
TEventList *targetCh; // for sorting events by their target positions
=======
// hold macroNo, targetTime, and targetPos for all target changer events
vector<tuple<unsigned int,double,unsigned int>> targetTimeList;
double targetTime; // this is the correct target time for each detector event
double prevTime = 0;
int prevTarget = 1;
double fullTime; // incorporates all 3 time stamps (fine, coarse, extended)
int targetPos = -1;


vector<bool> firstWaveform (8,false);

//TEventList *targetCh; // for sorting events by their target positions
>>>>>>> Stashed changes

bool text = false; // flag for producing text-file output apart from
                   // the default ROOT plots and tree filling
bool runlist = false; // flag indicating that runs should be read from
                      // runsToSort.txt, and NOT just the most recent file
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
    timetag *= 2; // timetag converted to ns from samples

    if(chNo==0 && firstWaveform[0])
    {
        // new macropulse
        macroNo++;
        //fill_n(evtNo.begin(),8,0);
    }
    timetagP[chNo] = timetag;

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
<<<<<<< Updated upstream
    out << "| run " << runNo << ", macro " << macroNo;
=======
    out << "| run " << runNo << ", macro " << left << setfill(' ') << setw(44) << "|" << endl;
    out << "| channel " << chNo;
>>>>>>> Stashed changes

    if(evtType==1)
    {
        out << left << setfill(' ') << setw(44) << ", DPP mode" << "|" << endl;
    }

    else if (evtType==2)
    {
        out << left << setfill(' ') << setw(44) << ", waveform mode " << "|" << endl;
    }

    else
    {
        out << left << setfill(' ') << setw(50) << ", ERROR DETERMINING MODE" << " |" << endl;
        cout << "Error: event type value out-of-range (DPP=1, waveform=2)" << endl;
        cout << "Event number = " << evtNo[chNo] << endl;
    }

    out << "| channel " << chNo << right << setfill(' ') << setw(52) << "|" << endl;

    temp << size << " bytes";
    out << "| size = " << left << setfill(' ') << setw(53) << temp.str() << "|" << endl;

    temp.str("");
    temp << 2*timetag << " ns";
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
<<<<<<< Updated upstream
    ev.runNo = std::stoi(runNo);
    ev.macroNo = macroNo;
    ev.evtNo = evtNo[chNo];
=======
    //ev.runNo = std::stoi(runNo);
    //ev.evtNo = evtNo[chNo];
>>>>>>> Stashed changes
    ev.chNo = chNo;
    ev.evtType = evtType;
    ev.extTime = extTime;
    ev.timetag = timetag;
    ev.fineTime = fineTime;
    ev.sgQ = sgQ;
    ev.lgQ = lgQ;
<<<<<<< Updated upstream
    ev.waveform = waveform;
    //ev.anProbe = anProbe;

    tree->Fill();

    // If an event comes from the target changer, we'll want to use it later
    // to group detector events based on target position. So segregate it into
    // its own TEventList that we'll loop over when we get detector events.
=======

    if (nE%10000 == 0)
    {
        ev.waveform = waveform;
    }
    //ev.anProbe = anProbe;

    tree->Fill();
>>>>>>> Stashed changes
}

void fillHistos()
{
<<<<<<< Updated upstream
    gDirectory->cd("/");
    gDirectory->GetDirectory(dirs[chNo].c_str())->cd();

    if (evtType == 1)
    {
        // DPP mode - fill DPP histograms

        outMacro = (TH1I*)(gDirectory->FindObject("outMacro"));
        outMacro->Fill(macroNo);

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
        outFT->Fill(fineTime);

        stringstream temp;
        temp << "evtNo " << evtNo[chNo];

        DPPWaveformsDir = (TDirectory*)gDirectory->FindObject("DPPWaveformsDir");
        DPPWaveformsDir->cd();

        if (evtNo[chNo]%100 == 0)
        {
            DPPWaveform = new TH1I(temp.str().c_str(),temp.str().c_str(),waveform.size(),0,waveform.size()*2);

            for(int i=0;i<waveform.size();i++)
            {
                DPPWaveform->SetBinContent(i,waveform[i]);
            }
        }

        fill_n(firstWaveform.begin(),8,true);
    }

    else if (evtType == 2)
    {
        // waveform mode - create full macropulse waveforms

        WaveWaveformsDir = (TDirectory*)gDirectory->FindObject("WaveWaveformsDir");
        WaveWaveformsDir->cd();

        switch (chNo)
        {
            case 0:

                if (firstWaveform[chNo])
                {
                    stringstream temp;
                    temp << "macropulse " << macroNo;

                    WaveWaveform0 = new TH1I(temp.str().c_str(),temp.str().c_str(),(WAVELET_NO+1)*waveform.size(),0,(WAVELET_NO+1)*waveform.size()*2);

                    waveformStart[chNo] = timetag;

                    firstWaveform[chNo] = false;
                }

                for(int i=0;i<waveform.size();i++)
                {
                    WaveWaveform0->SetBinContent(i+timetag-waveformStart[chNo],waveform[i]);
                }

                break;

            case 2:

                if (firstWaveform[chNo])
                {
                    stringstream temp;
                    temp << "macropulse " << macroNo;

                    WaveWaveform2 = new TH1I(temp.str().c_str(),temp.str().c_str(),(WAVELET_NO+1)*waveform.size(),0,(WAVELET_NO+1)*waveform.size()*2);

                    waveformStart[chNo] = timetag;

                    firstWaveform[chNo] = false;
                }

                for(int i=0;i<waveform.size();i++)
                {
                    WaveWaveform2->SetBinContent(i+timetag-waveformStart[chNo],waveform[i]);
                }

                break;

            case 4:

                if (firstWaveform[chNo])
                {
                    stringstream temp;
                    temp << "macropulse " << macroNo;

                    WaveWaveform4 = new TH1I(temp.str().c_str(),temp.str().c_str(),(WAVELET_NO+1)*waveform.size(),0,(WAVELET_NO+1)*waveform.size()*2);

                    waveformStart[chNo] = timetag;

                    firstWaveform[chNo] = false;
                }

                for(int i=0;i<waveform.size();i++)
                {
                    WaveWaveform4->SetBinContent(i+timetag-waveformStart[chNo],waveform[i]);
                }

                break;

            case 6:

                if (firstWaveform[chNo])
                {
                    stringstream temp;
                    temp << "macropulse " << macroNo;

                    WaveWaveform6 = new TH1I(temp.str().c_str(),temp.str().c_str(),(WAVELET_NO+1)*waveform.size(),0,(WAVELET_NO+1)*waveform.size()*2);

                    waveformStart[chNo] = timetag;

                    firstWaveform[chNo] = false;
                }

                for(int i=0;i<waveform.size();i++)
                {
                    WaveWaveform6->SetBinContent(i+timetag-waveformStart[chNo],waveform[i]);
                }

                break;

            case 7:

                if (firstWaveform[chNo])
                {
                    stringstream temp;
                    temp << "macropulse " << macroNo;

                    WaveWaveform7 = new TH1I(temp.str().c_str(),temp.str().c_str(),(WAVELET_NO+1)*waveform.size(),0,(WAVELET_NO+1)*waveform.size()*2);

                    waveformStart[chNo] = timetag;

                    firstWaveform[chNo] = false;
                }

                for(int i=0;i<waveform.size();i++)
                {
                    WaveWaveform7->SetBinContent(i+timetag-waveformStart[chNo],waveform[i]);
=======
    tree->SetBranchAddress("fineTime",&fineTime);
    tree->SetBranchAddress("waveform",&waveform);

    int totalEntries = tree->GetEntries();

    double gTimes;
    int gCounter;
    TH1I* gammas = new TH1I("gammas","gammas",100,70,110);

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

    tree->Draw(">>targetChEvents","chNo==0","entrylist");
    tree->Draw(">>monitorEvents","chNo==2","entrylist");
    tree->Draw(">>detSEvents","chNo==4","entrylist");
    tree->Draw(">>detLEvents","chNo==6","entrylist");
    tree->Draw(">>detREvents","chNo==7","entrylist");

    TEntryList* targetChEvents = (TEntryList*)gDirectory->Get("targetChEvents");
    TEntryList* monitorEvents = (TEntryList*)gDirectory->Get("monitorEvents");
    TEntryList* detSEvents = (TEntryList*)gDirectory->Get("detSEvents");
    TEntryList* detLEvents = (TEntryList*)gDirectory->Get("detLEvents");
    TEntryList* detREvents = (TEntryList*)gDirectory->Get("detREvents");

    vector<TEntryList*> channelList;
    channelList.push_back(targetChEvents);
    channelList.push_back(monitorEvents);
    channelList.push_back(detSEvents);
    channelList.push_back(detLEvents);
    channelList.push_back(detREvents);

    for(int j = 0; j<channelList.size(); j++)
    {
        // loop through tree once per channel number 
        cout << "Populating " << dirs[j] << " histograms..." << endl;
        nE = 0;

        gDirectory->cd("/");
        gDirectory->GetDirectory(dirs[j].c_str())->cd();

        DPPWaveformsDir = (TDirectory*)gDirectory->FindObject("DPPWaveformsDir");
        WaveWaveformsDir = (TDirectory*)gDirectory->FindObject("WaveWaveformsDir");

        cout << "found dpp/waveforms directories" << endl;

        int targetCounter = 0;
        targetTime = get<1>(targetTimeList[targetCounter]);

        tree->SetEntryList(channelList[j]);
        totalEntries = tree->GetEntries();

        for (int i=0; i<totalEntries; i++)

        {
            tree->GetEntry(i);

        
            if (chNo == chNoList[j])
            {
                if (evtType == 1)
                {

                    // prepare for filling basic histos
                    fullTime = timetag+pow(2,32)*extTime;

                    if (chNo == 4 || chNo == 6 || chNo == 7)
                    {
                        fullTime += fineTime*(2./1024.);
                    }

                    while (fullTime-targetTime+TIME_OFFSET > 650000)
                    {
                        // if it's been too long since the last target changer event,
                        // step to the next target changer event
                        targetCounter++;
                        targetTime = get<1>(targetTimeList[targetCounter]);
                        targetPos = get<2>(targetTimeList[targetCounter]);
                        fill(evtNo.begin(),evtNo.end(),0);
                    }

                    // if event has associated target changer event, fill DPP histo
                    if (fullTime-targetTime+TIME_OFFSET >= 0) 
                    {
                        outMacro = (TH1I*)(gDirectory->FindObject("outMacro"));
                        outMacro->Fill(get<0>(targetTimeList[targetCounter]));

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
                        outFT->Fill(fineTime + 16*rng->Rndm());

                        if (nE%10000 == 0)
                        {
                            DPPWaveformsDir->cd();

                            stringstream temp;
                            temp << "macroNo " << get<0>(targetTimeList[targetCounter]) << "evtNo " << evtNo[chNo];
                            DPPWaveform = new TH1I(temp.str().c_str(),temp.str().c_str(),waveform.size(),0,waveform.size()*2);

                            for(int i=0;i<waveform.size();i++)
                            {
                                DPPWaveform->SetBinContent(i,waveform[i]);
                            }
                        }
                    }

                    fill(waveformStart.begin(),waveformStart.end(),-1000000); // prepare for next waveform mode

                    if (chNo==2 || chNo==4 || chNo==6 || chNo==7)
                    {
                        float trueTime = fmod(((double)timetag+((double)fineTime*2./1024.)-targetTime-TIME_OFFSET),MICRO_PERIOD);

                        // convert trueTime into neutron velocity based on flight path distance
                        float velocity = pow(10,7)*FLIGHT_DISTANCE/trueTime; // in meters/sec 

                        // convert velocity to relativistic kinetic energy
                        float rKE = (pow((1-pow((velocity/C),2)),-0.5)-1)*NEUTRON_MASS; // in MeV

                        if (trueTime > 100 || (chNo==2 && trueTime > 50)) // gate disallowing gammas
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

                        if (targetPos > 0)
                        {
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
                        } 
                    }
                }

                else if (evtType == 2)
                {
                    WaveWaveformsDir->cd();

                    TH1I* waveformHolder;

                    if (fullTime > waveformStart[chNo] + 1000000)
                    {
                        // new macropulse in waveform mode - create new plot

                        stringstream temp;
                        temp << "macropulse " << get<0>(targetTimeList[targetCounter])+1;
                        waveformStart[chNo] = fullTime;

                        switch (chNo)
                        {
                            case 0:
                                WaveWaveform0 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);
                                waveformHolder = WaveWaveform0;

                                targetCounter++; // increment target changer number in list to prepare for next DPP mode

                                break;

                            case 2:
                                WaveWaveform2 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);
                                waveformHolder = WaveWaveform2;
                                break;

                            case 4:
                                WaveWaveform4 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);
                                waveformHolder = WaveWaveform4;
                                break;

                            case 6:
                                WaveWaveform6 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);
                                waveformHolder = WaveWaveform6;
                                break;

                            case 7:
                                WaveWaveform7 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);
                                waveformHolder = WaveWaveform7;
                                break;
                        }
                    }

                    for(int i=0;i<waveform.size();i++)
                    {
                        waveformHolder->SetBinContent(i+(fullTime-waveformStart[chNo])/2,waveform[i]);
                    }
                }

                nE++;
                evtNo[chNo]++;

                if(nE%100==0)
                {
                    cout << nE << " events\r";
                    fflush(stdout);
>>>>>>> Stashed changes
                }

                break;
        }
<<<<<<< Updated upstream
=======
        cout << endl;
>>>>>>> Stashed changes
    }
}

void calculateCS()
{
    // create cross-section histograms

<<<<<<< Updated upstream
    gDirectory->cd("/");
=======
    //fill(targetTimeList.begin(),targetTimeList.end(),make_pair(0,0)); // clear the target time list

>>>>>>> Stashed changes

    TH1I *blank = new TH1I("blank","blank",1000,0,800);
    TH1I *carbonS = new TH1I("carbonS","carbonS",1000,0,800);
    TH1I *carbonL = new TH1I("carbonL","carbonL",1000,0,800);
    TH1I *Sn112 = new TH1I("Sn112","Sn112",1000,0,800);
    TH1I *NatSn = new TH1I("NatSn","NatSn",1000,0,800);
    TH1I *Sn124 = new TH1I("Sn124","Sn124",1000,0,800);

    TH1I *TOF = new TH1I("TOF","Time of flight",1000,0,2000);

    // link tree branches to variables-to-read
<<<<<<< Updated upstream
    tree->SetBranchAddress("macroNo",&macroNo);
=======
>>>>>>> Stashed changes
    tree->SetBranchAddress("chNo",&chNo);
    tree->SetBranchAddress("lgQ",&lgQ);
    tree->SetBranchAddress("timetag",&timetag);
    tree->SetBranchAddress("fineTime",&fineTime);
    tree->SetBranchAddress("extTime",&extTime);

<<<<<<< Updated upstream
    // create a list of all target-changer events for sub-looping
    tree->Draw(">>targetCh","chNo == 0","entrylist");
    TEntryList *targetCh;
    gDirectory->GetObject("targetCh",targetCh);
    tree->SetEntryList(targetCh);
    int targetEntries = tree->GetEntries();

    int targetPos = -1;

    //int previousTargetPos = -1;
=======
    int totalEntries = tree->GetEntries();
>>>>>>> Stashed changes
    
    //map<long,int> targetMap;

    /*for (int i=0; i<targetEntries; i++)
    {
        tree->GetEntry(i);
        
        // assign an integral target position based on lgQ of target changer signal
        if (lgQ>tarGate[0] && lgQ<tarGate[1])
        {
            currentTargetPos = 1;
        }

        else if (lgQ>tarGate[2] && lgQ<tarGate[3])
        {
<<<<<<< Updated upstream
            currentTargetPos = 2;
        }

        else if (lgQ>tarGate[4] && lgQ<tarGate[5])
        {
            currentTargetPos = 3;
        }

        else if (lgQ>tarGate[6] && lgQ<tarGate[7])
        {
            currentTargetPos = 4;
        }

        else if (lgQ>tarGate[8] && lgQ<tarGate[9])
        {
            currentTargetPos = 5;
        }

        else if (lgQ>tarGate[10] && lgQ<tarGate[11])
        {
            currentTargetPos = 6;
        }

        if (currentTargetPos != previousTargetPos)
        {
            // target position change detected

            long targetTime = pow(2,32)*extTime + timetag + fineTime*(2000/1024); // in ns
            targetMap.insert(pair<long,int>(targetTime,currentTargetPos)); // gives time of target change and new position
            previousTargetPos = currentTargetPos;
        }
    }*/

    // create a list of all events in the tree
    int totalEntries = tree->GetEntries();
    tree->Draw(">>total","","entrylist");
    TEntryList *total;
    gDirectory->GetObject("total",total);
        
    tree->SetEntryList(total);

    for (int i=0; i<totalEntries/10; i++)
    {
        tree->GetEntry(i);

        if (chNo == 4 || chNo == 6 || chNo == 7) // detector event; continue
        {
            int detMacro = macroNo; // save the detector event's data for the energy calculation
            int detExtTime = extTime;
            int detCoarseTime = timetag;
            int detFineTime = fineTime;

            for (int j=i; j>0; j--)
=======
            if (chNo == 2 || chNo == 4 || chNo == 6 || chNo == 7)
>>>>>>> Stashed changes
            {
                tree->GetEntry(j);

<<<<<<< Updated upstream
                //cout << "detMacro is " << detMacro << ", targetMacro is " << macroNo << endl;

                if (chNo == 0 && detMacro == macroNo) // same macropulse; found target changer position
                {
                    // assign an integral target position based on lgQ of target changer signal
                    if (lgQ>tarGate[0] && lgQ<tarGate[1])
                    {
                        targetPos = 1;
=======
                for(vector<tuple<unsigned int,double,unsigned int>>::reverse_iterator it = targetTimeList.rbegin(); it != targetTimeList.rend(); ++it)
                {
                    if ((timetag+pow(2,32)*extTime) > get<1>(*it) && (timetag+pow(2,32)*extTime) <get<1>(*it)+650000) 
                    {
                        // found a previous target time within 1 ms 
                        targetTime = get<1>(*it);
                        targetPos = get<2>(*it);
                        break;
>>>>>>> Stashed changes
                    }

                    if (lgQ>tarGate[2] && lgQ<tarGate[3])
                    {
                        targetPos = 2;
                    }

                    if (lgQ>tarGate[4] && lgQ<tarGate[5])
                    {
                        targetPos = 3;
                    }

                    if (lgQ>tarGate[6] && lgQ<tarGate[7])
                    {
                        targetPos = 4;
                    }

                    if (lgQ>tarGate[8] && lgQ<tarGate[9])
                    {
                        targetPos = 5;
                    }

                    if (lgQ>tarGate[10] && lgQ<tarGate[11])
                    {
                        targetPos = 6;
                    }

                    // calculate a targetTime-calibrated time-of-flight of neutrons
                    // trueTime is neutron time-of-flight since micropulse start
                    float trueTime = fmod((detCoarseTime+(detFineTime*2/1024)-timetag+TIME_OFFSET),MICRO_PERIOD);

                    float fakeTime = fmod((detCoarseTime-timetag),MICRO_PERIOD);

                    // convert trueTime into neutron velocity based on flight path distance
                    float velocity = pow(10,9)*FLIGHT_DISTANCE/trueTime; // in meters/sec 

                    // convert velocity to relativistic kinetic energy
                    float rKE = (pow((1-pow((velocity/C),2)),-0.5)-1)*NEUTRON_MASS*pow(C,2); // in MeV

                    switch (targetPos)
                    {
                        case 1:
                            // BLANK
                            blank->Fill(rKE);
                            break;
                        case 2:
                            // SHORT CARBON
                            carbonS->Fill(rKE);
                            break;
                        case 3:
                            // LONG CARBON
                            carbonL->Fill(rKE);
                            break;
                        case 4:
                            // Sn112
                            Sn112->Fill(rKE);
                            break;
                        case 5:
                            // Natural Sn
                            NatSn->Fill(rKE);
                            break;
                        case 6:
                            // Sn124
                            Sn124->Fill(rKE);
                            break;
                    }

<<<<<<< Updated upstream
                    TOF->Fill(fakeTime);

                    if (i%100 == 0)
                    {
                        cout << "Populated " << i << " events in cross-section histograms\r";
                        fflush(stdout);
                    }
=======
                if (targetPos > 0)
                {
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
                }
>>>>>>> Stashed changes

                    break;
                }
            }
        }
    }

    cout << endl << endl;
}


void processRun(string evtname)
{
    // attempt to process the event file
    ifstream evtfile;
    evtfile.open(evtname,ios::binary);

    // for text file displaying events, organized by input channel to make sense of the data
    ofstream totalOut;
    ofstream targetChangerOut;
    ofstream monitorOut;
    ofstream detectorLOut;
    ofstream detectorROut;
    ofstream detectorTOut;

    if(text)
    {
        totalOut.open("textSort/allChannels.txt");
        targetChangerOut.open("textSort/targetChanger.txt");
        monitorOut.open("textSort/monitor.txt");
        detectorLOut.open("textSort/detectorL.txt");
        detectorROut.open("textSort/detectorR.txt");
        detectorTOut.open("textSort/detectorT.txt");
    }

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

        // start looping through the evtfile for events
<<<<<<< Updated upstream
        while(!evtfile.eof())
=======

        int looplimit = 0;

        while(!evtfile.eof() && looplimit<1000000)
>>>>>>> Stashed changes
        {
            // get channel number from the event header
            chNo = readHeader(evtfile);

            readBody(evtfile); // extract data from the event body
            fillTree(); // fill the tree with event data
            fillHistos(); // fill the histos with event data

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

<<<<<<< Updated upstream
=======
void populateMacros()
{
    cout << "Creating macropulse number and time list" << endl;

    // link tree branches to variables-to-read
    tree->SetBranchAddress("chNo",&chNo);
    tree->SetBranchAddress("lgQ",&lgQ);
    tree->SetBranchAddress("timetag",&timetag);
    tree->SetBranchAddress("extTime",&extTime);

    int totalEntries = tree->GetEntries();

    unsigned int macroNo = 0; // keep track of target changer macropulse number

    for (int i=0; i<totalEntries; i++)
    {
        tree->GetEntry(i);

        // Refresh targetCh time
        if (chNo==0) 
        {
            if (evtType==1)
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

                // new macropulse in DPP mode
                macroNo++;
                fill_n(evtNo.begin(),8,0); // clear the evtNo counter for each channel

                targetTimeList.push_back(make_tuple(macroNo,timetag+pow(2,32)*extTime,targetPos));
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

            prevTarget = evtType;

            if (macroNo%100==0)
            {
                cout << "Macro number " << macroNo << ", fullTime " << timetag+pow(2,32)*extTime << "\r";
                fflush(stdout);
            }

        }
    }
    cout << endl;
}

>>>>>>> Stashed changes
int main(int argc, char* argv[])
{

    /*************************************************************************/
    /* read flags to set the sorting mode                                    */
    /*************************************************************************/ 

    if (argc > 1) // flags detected
    {
        if (std::string(argv[2]) == "true")
        {
            // produce text files for each channel containing all event
            // data from the input file. This will significantly increase
            // processing time for this program
            text = true;
        }

        if (string(argv[3]) == "true")
        {
            // produce cross-section plots based on target position
            cs = true;
        }
    }

    runNo = argv[1];

    // Create a tree for this run
    TFile *file;

    stringstream treeName;
    stringstream fileName;
<<<<<<< Updated upstream
    treeName << "run" << runNo; 
    fileName << treeName.str() << ".root";

    file = new TFile(fileName.str().c_str(),"RECREATE");

    tree = new TTree("tree","");
    cout << "Created ROOT tree " << treeName.str() << endl;

    tree->Branch("runNo",&ev.runNo,"runNo/i");
    tree->Branch("macroNo",&ev.macroNo,"macroNo/i");
    tree->Branch("evtNo",&ev.evtNo,"evtNo/i");
    tree->Branch("chNo",&ev.chNo,"chNo/i");
    tree->Branch("evtType",&ev.evtType,"evtType/i");
    tree->Branch("extTime",&ev.extTime,"extTime/i");
    tree->Branch("timetag",&ev.timetag,"timetag/i");
    tree->Branch("fineTime",&ev.fineTime,"fineTime/i");
    tree->Branch("sgQ",&ev.sgQ,"sgQ/i");
    tree->Branch("lgQ",&ev.lgQ,"lgQ/i");
    tree->Branch("waveform",&ev.waveform);

    //ENABLE for diagnostic analog probe
    //tree->Branch("anProbe",&ev.anProbe);
        
    for(int i=0; i<8; i++)
=======
    treeName << runDir << "-" << runNo; 
    fileName << "/media/ExternalDrive1/analysis/" << runDir << "/" << treeName.str() << ".root";

    file = new TFile(fileName.str().c_str(),"RECREATE");

    if(file->Get("tree"))
>>>>>>> Stashed changes
    {
        vector<TH1I*> tempVec;
        histos.push_back(tempVec); // create sub-vector for this channel

<<<<<<< Updated upstream
        if (dirs[i].compare("") != 0) // valid channel
        {
            gDirectory->mkdir(dirs[i].c_str(),dirs[i].c_str());
            gDirectory->GetDirectory(dirs[i].c_str())->cd();

            // instantiate histograms
=======
    else
    {
        tree = new TTree("tree","");
        cout << "Created ROOT tree " << treeName.str() << endl;

        //tree->Branch("runNo",&ev.runNo,"runNo/i");
        //tree->Branch("evtNo",&ev.evtNo,"evtNo/i");
        tree->Branch("chNo",&ev.chNo,"chNo/i");
        tree->Branch("evtType",&ev.evtType,"evtType/i");
        tree->Branch("extTime",&ev.extTime,"extTime/i");
        tree->Branch("timetag",&ev.timetag,"timetag/i");
        tree->Branch("fineTime",&ev.fineTime,"fineTime/i");
        tree->Branch("sgQ",&ev.sgQ,"sgQ/i");
        tree->Branch("lgQ",&ev.lgQ,"lgQ/i");
        tree->Branch("waveform",&ev.waveform);

        //ENABLE for diagnostic analog probe
        //tree->Branch("anProbe",&ev.anProbe);

        for(int i=0; i<5; i++)
        {
            vector<TH1I*> tempVec;
            histos.push_back(tempVec); // create sub-vector for this channel

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
>>>>>>> Stashed changes

            histos.back().push_back(new TH1I("outMacro","outMacro",100000,0,10000000));
            histos.back().push_back(new TH1I("outEvt","outEvt",500,0,1000));
            histos.back().push_back(new TH1I("outExtTime","outExtTime",1000,0,1000));
            histos.back().push_back(new TH1I("outTime","outTime",2500000,0,2500000000));
            histos.back().push_back(new TH1I("outSGQ","outSGQ",1024,0,70000));
            histos.back().push_back(new TH1I("outLGQ","outLGQ",1024,0,70000));
            histos.back().push_back(new TH1I("outFT","outFT",1023,0,1023));

            gDirectory->mkdir("DPPWaveformsDir","raw DPP waveforms");
            gDirectory->mkdir("WaveWaveformsDir","concatenated waveform waveforms");

            gDirectory->cd("/");
        }
    }

<<<<<<< Updated upstream
    stringstream runName;
    runName << "../output/run" << runNo << ".evt";
    processRun(runName.str());

    if (cs)
    {
        calculateCS(); // produce histograms showing cross-sections for each target
    }

=======
>>>>>>> Stashed changes
    file->Write();

}
