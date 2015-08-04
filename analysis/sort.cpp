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
unsigned long timetag, timetagP = 0;

unsigned int Ne = 0; // counter for the total number of events in the input file
unsigned int Nwavelets = 0; // counter for the number of MIXED mode wavelets in the input file 
unsigned int Nwaveforms = 0; // counter for the number of waveform events in the input file 

// Create a ROOT tree for event data
TTree* tree;
TH1S* outDPP;

vector<TH1S*> listWaveforms;
vector<TH1S*> listWavelets; 

/*struct DPPevent {
    unsigned int timetag, sgQ, lgQ, baseline;
} de;

tree->Branch("event",&de.timetag,"time/I:sgQ:lgQ:baseline");

outDPP = new TH1S("outDPP","outDPP",1024,0,60000);
*/

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

// prints EVENT HEADER to the appropriate text file
// 'out' points to a channel-specific text file
void printHeader(ofstream& out)
{

    // use for formatting strings with units
    stringstream temp;

    out << setfill('*') << setw(63) << "*" << endl;
    out << "| EVENT " << left << setfill(' ') << setw(54) << Ne << "|" << endl;
    out << "|" << right << setfill('-') << setw(62) << "|" << endl;
    out << "| channel #" << channel;

    if(evtype==1)
    {
        out << left << setfill(' ') << setw(50) << ", DPP mode" << "|" << endl;
    }
    else if (evtype==2)
    {
        out << left << setfill(' ') << setw(50) << ", waveform mode " << "|" << endl;
    }
    else
    {
        out << left << setfill(' ') << setw(50) << ", ERROR DETERMINING MODE" << " |" << endl;
        cout << "Error: event type value out-of-range (DPP=1, waveform=2)" << endl;
        cout << "Event number = " << Ne << endl;
    }

    temp << size << " bytes";
    out << "| size = " << left << setfill(' ') << setw(53) << temp.str() << "|" << endl;

    temp.str("");
    temp << timetag << " ns";
    out << "| timetag = " << left << setw(50) << temp.str() << "|" << endl;

    out << "|" << right << setfill('-') << setw(62) << "|" << endl;

}

// unpacks EVENT BODY data into ROOT histograms and text output
// 'out' points to a channel-specific text file
void unpack(ifstream& evtfile, ofstream& out)
{
    ostringstream histName;
    
    // For DPP EVENT BODY unpacking 
    if(evtype==1)
    {

        // extras is a 32-bit word whose content varies with the value of EXTRA_SELECT.
        // [DEFAULT]    0: extended timestamp (bits 16-31) and baseline*4 (bits 0-15)
        //              1: extended timestamp (bits 16-31) and flags (bits 0-15?)
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
                // flag documentation
                out << "| flags = " << left << setfill(' ') << setw(49) << extras1 << "|" << endl;
                break; 

            case 2:
                temp << extras2 << " *8.59 s";
                out << "| extended time stamp = " << left << setfill(' ') << setw(34) << temp.str() << "|" << endl;
                // flag documentation
                out << "| flags = " << left << setfill(' ') << setw(49) << extras1 << "|" << endl;
                // fine time stamp documentation
                out << "| fine time stamp = " << left << setfill(' ') << setw(49) << extras1 << "|" << endl;
                break;

            case 3:
                // pulse peak value documentation
                out << "| pulse peak value" << left << setfill(' ') << setw(49) << extras1 << "|" << endl;
                break;

            case 5:
                // PZC and NZC
                out << "| PZC" << left << setfill(' ') << setw(49) << extras2 << "|" << endl;
                out << "| NZC" << left << setfill(' ') << setw(49) << extras1 << "|" << endl;
                break;

            case 7:
                // fixed value of 0x12345678
                const unsigned int extras = (extras2 << 16) | extras1;
                out << "| " << left << setfill(' ') << setw(49) << extras << "|" << endl;
                break;
        }

        // sgQ is the short gate integrated charge, in digitizer units 
        unsigned short sgQ = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        out << "| short gate charge = " << left << setw(40) << sgQ << "|" << endl;

        // lgQ is the long gate integrated charge, in digitizer units 
        unsigned short lgQ = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        //outDPP->Fill(lgQ);
        out << "| long gate charge = " << left << setw(41) << lgQ << "|" << endl;

        // baseline is the baseline level, frozen at the trigger time
/*        unsigned short baseline = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        out << "| baseline level = " << left << setw(42) << baseline << " |" << endl;
*/
        // puRej is a pile-up rejection flag (1 if pile-up detected? 0 if not?)
        unsigned short puRej = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        //out << "| pile-up detected = " << left << setw(40) << puRej << " |" << endl;

        // nSamp is the number of waveform samples that follow (in LIST mode, this is 0)
        unsigned short nSamp1 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        unsigned short nSamp2 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        const unsigned int nSamp = (nSamp2 << 16) | nSamp1;
        temp.str("");
        temp << nSamp << " samples";
        out << "| waveform length = " << left << setw(41) << temp.str() << " |" << endl;
        out << "|" << right << setfill('-') << setw(62) << "|" << endl;

        if(nSamp > 0)
            // We must be in MIXED mode; time to output wavelet.
        {
            histName.str("");
            histName << "outWavelet" << Ne;
            string tempHist = histName.str();

            listWavelets.push_back(new TH1S(tempHist.c_str(),tempHist.c_str(),nSamp,0,nSamp*2));

            out << left << setfill(' ') << setw(62) << "| Waveform samples" << "|" << endl;
            out << "|" << right << setfill(' ') << setw(62) << "|" << endl;

            for(int i=0;i<nSamp;i++)
            {
                listWavelets[Nwavelets]->SetBinContent(i,buffer[0]);

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
                evtfile.read((char*)buffer,BufferBytes);

                if(i%10==9 || i==nSamp-1)
                {
                    out << left << setw(62) << temp.str() << "|" << endl;
                } 

            }
            
            // done with this event; increment the wavelet counter to get ready for the next
            // wavelet.
            Nwavelets++;
        }

        // Fill the root tree with data extracted from the event for later analysis
        /*de.timetag = timetag;
        de.sgQ = sgQ;
        de.lgQ = lgQ;
        de.baseline = baseline;

        tree->Fill();
        */
    }

    else if(evtype==2)
    {

        unsigned short nSamp1 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        unsigned short nSamp2 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        unsigned long nSamp = (nSamp2 << 16) | nSamp1;

        stringstream temp;
        temp << nSamp << " samples";
        out << "| waveform length = " << left << setfill(' ') << setw(41) << temp.str() << " |" << endl;
        out << "|" << right << setfill('-') << setw(62) << "|" << endl;

        histName.str("");
        histName << "outWaveform" << Ne;
        string tempHist = histName.str();

        listWaveforms.push_back(new TH1S(tempHist.c_str(),tempHist.c_str(),nSamp,0,nSamp*2));

        for(int i=0;i<nSamp;i++)
        {

            if(i%1000 == 0)
            {
                out << "|" << right << setfill(' ') << setw(62) << "|" << endl;
                stringstream temp;
                temp << "Samples " << i << "-" << i+1000;
                out << "| " << left << setw(60) << temp.str() << "|" << endl;
}

            if(i%100 == 0)
            {
                out << "|" << right << setw(62) << "|" << endl;
            }

            if(i%10 == 0)
            {
                out << "|";
            } 

            listWaveforms[Nwaveforms]->SetBinContent(i,buffer[0]);
            out << right << setfill(' ') << setw(6) << buffer[0];
            evtfile.read((char*)buffer,BufferBytes);

            if(i%10 == 9)
            {
                out << " |" << endl;
            } 
        }

        // done with this event; increment the wavelet counter to get ready for the next
        // wavelet.

        out << "\n";
        Nwaveforms++;
    }

    else
    {
        cout << "ERROR: unknown value for event type" << endl;
        return;
    }

    out << setfill('*') << setw(63) << "*" << endl;
    out << endl;

}

int main()
{

    // create a ROOT file for holding events, with one directory per input channel
    TFile *file; 
    file = new TFile("events.root","RECREATE");
    file->cd();
    
    // create a ROOT tree for holding event properties
    tree = new TTree("tree","");
    tree->SetAutoSave(0);

    TDirectoryFile *targetChangerDir, *monitorDir, *detLDir, *detRDir, *detTDir;
    targetChangerDir = new TDirectoryFile("targetChanger","Target Changer");
    monitorDir = new TDirectoryFile("monitor","Monitor");
    detLDir = new TDirectoryFile("detL","Detector left)");
    detRDir = new TDirectoryFile("detR","Detector right");
    detTDir = new TDirectoryFile("detT","Detector sum");
    // TDirectory *dirSub; to be uncommented and used for sub-directories, if desired
 
    // create text output for holding events, with file per input channel
    ofstream totalOut ("sorted/allChannels.txt");
    ofstream targetChangerOut ("sorted/targetChanger.txt");
    ofstream monitorOut ("sorted/monitor.txt");
    ofstream detectorLOut ("sorted/detectorL.txt");
    ofstream detectorROut ("sorted/detectorR.txt");
    ofstream detectorTOut ("sorted/detectorT.txt");

    // create text output to examine time differences from one event to the next
    ofstream timeDiff ("sorted/timeDiff.txt");

    // open the event file
    ifstream evtfile;
    string name = "../output/wutest.evt";

    evtfile.clear();
    evtfile.open(name.c_str(),ios::binary);      

    if (evtfile.bad()) cout << "bad" << endl;
    if (evtfile.is_open()) cout << "open" << endl;
    if (evtfile.good()) cout << "good" << endl;
    if (evtfile.fail()) cout << "fail" << endl;

    evtfile.read((char*)buffer,BufferBytes);

    point = buffer;

    // start looping through the evtfile for events
    while(!evtfile.eof())
    {
        // get channel number from the event header
        int channelNum = readHeader(evtfile);

        // print the header to an allChannels.txt
        printHeader(totalOut);

        // print the time difference between adjacent events to timeDiff.txt
        if(timetag>timetagP)
        {
            timeDiff << timetag - timetagP << endl;
        }
        timetagP = timetag;

        // funnel event data differently based on channel
        switch (channelNum)
        {
            case 0: 
                // target-changer data
                targetChangerDir->cd();
                printHeader(targetChangerOut);
                unpack(evtfile, targetChangerOut);
                break;

            case 1:
                break;

            case 2:
                // monitor data
                monitorDir->cd();
                printHeader(monitorOut);
                unpack(evtfile, monitorOut);
                break;

            case 3:
                break;

            case 4:
                // Detector data (assume detectors T'd together)
                detTDir->cd();
                printHeader(detectorTOut);
                unpack(evtfile, detectorTOut);
                break;

            case 5:
                break;

            case 6:
                // Detector data (left detector only)
                detLDir->cd();
                printHeader(detectorLOut);
                unpack(evtfile, detectorLOut);
                break;

            case 7:
                // Detector data (right detector only)
                detRDir->cd();
                printHeader(detectorROut);
                unpack(evtfile, detectorROut);
                break;

            default:
                cout << "ERROR: unknown value for channel type" << endl;
                return 0;
                break;
        }

        // Event finished
        Ne++;
    }
    
    // Input file finished
    cout << "Finished processing event file" << endl;
    cout << "Total events: " << Ne << endl;

    file->Write();
    tree->Write();
}


/*struct MIXEDevent {
  vector<unsigned int> wavelet;
  };

  MIXEDevent me;

  for(int i=0;i<nSamp;i++)
  {
  me.wavelet.push_back(buffer[0]);
  evtfile.read((char*)buffer,BufferBytes);
  }


  }*/
