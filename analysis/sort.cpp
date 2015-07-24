// C++ file to read CAEN event files 

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TH1S.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

// read the header for the event to determine how it should be processed
// returns the channel number of the event to the main processing function

int const BufferWords = 1;
int const BufferBytes = BufferWords*2;
unsigned short buffer[BufferWords];
unsigned short *point;
unsigned long size;
unsigned long evtype;
unsigned long channel;
unsigned long timetag;

unsigned int Ne = 0; // counter for the number of events in the buffer
unsigned int Nwavelets = 0; // counter for the number of MIXED mode wavelets in the buffer
unsigned int Nwaveforms = 0; // counter for the number of waveform events in the buffer


ofstream timeOut ("timeTags.txt");
ofstream totalOut ("output.txt");
ofstream targetChangerOut ("targetChangerOut.txt");
ofstream monitorOut ("monitorOut.txt");
ofstream detectorOut ("detectorOut.txt");

ostringstream histName;
ifstream evtfile;
string name = "../output/wutest.evt";

vector<TH1S*> listWaveforms;
vector<TH1S*> listWavelets; 

struct DPPevent {
    unsigned int timetag, sgQ, lgQ, baseline;
};

DPPevent de;

TFile* f; 
TTree* tree;

TH1S* outDPP;

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

// prints the data from the event header to the appropriate outfile
void headerPrint(ofstream& out)
{
    out << "size = " << size << endl;

    if(evtype==1)
    {
        out << "event type = DPP" << endl;
    }
    else if (evtype==2)
    {
        out << "event type = waveform" << endl;
    }
    else
    {
        cout << "Error: event type value out-of-range (DPP=1, waveform=2)" << endl;
    }

    out << "channel # = " << channel << endl;
    out << "timetag = " << timetag << endl;

}

// unpack takes data from an output.evt into readable ROOT or .txt output files, based on channel
// the parameter indicates how the data should be interpreted, i.e., "detector" or "target changer"

void unpack(ifstream& evtfile, ofstream& out)
{
    if(evtype==1)
        // DPP mode detected
    {

        // zcTime is the interpolated zero-crossing time, in picoseconds
        unsigned short zcTime = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        out << "zero crossing Time = " << zcTime << endl;

        // sgQ is the short gate integrated charge, in digitizer units 
        unsigned short sgQ = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        out << "short gate charge = " << sgQ << endl;

        // lgQ is the long gate integrated charge, in digitizer units 
        unsigned short lgQ = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        outDPP->Fill(lgQ);
        out << "long gate charge = " << lgQ << endl;

        // baseline is the baseline level, frozen at the trigger time
        unsigned short baseline = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        out << "baseline level = " << baseline << endl;

        // puRej is a pile-up rejection flag (1 if pile-up detected? 0 if not?)
        unsigned short puRej = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        out << "pile-up detected = " << puRej << endl;

        // nSamp is the number of waveform samples that follow (in LIST mode, this is 0)
        unsigned short nSamp1 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        unsigned short nSamp2 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        const unsigned int nSamp = (nSamp2 << 16) | nSamp1;
        out << "number of waveform samples = " << nSamp << endl;

        if(nSamp > 0)
            // We must be in MIXED mode; time to output wavelet.
        {
            out << "MIXED event wavelet data";

            histName.str("");
            histName << "outWavelet" << Nwavelets;
            string tempHist = histName.str();

            listWavelets.push_back(new TH1S(tempHist.c_str(),tempHist.c_str(),nSamp,0,nSamp*2));

            for(int i=0;i<nSamp;i++)
            {
                if(i%8 == 0)
                {
                    totalOut << "\n";
                }

                listWavelets[Nwavelets]->SetBinContent(i,buffer[0]);
                out << buffer[0] << " ";
                evtfile.read((char*)buffer,BufferBytes);

            }
            // done with this event; increment the wavelet counter to get ready for the next
            // wavelet.
            out << "\n\n";
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
        out << "Waveform event data" << endl;

        unsigned short nSamples1 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        unsigned short nSamples2 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        unsigned long nSamples = (nSamples2 << 16) | nSamples1;
        out
            << "nSamples = " << nSamples << endl;

        //unsigned short waveform[nSamples];

        histName.str("");
        histName << "outWaveform" << Nwaveforms;
        string tempHist = histName.str();

        listWaveforms.push_back(new TH1S(tempHist.c_str(),tempHist.c_str(),nSamples,0,nSamples*2));

        for(int i=0;i<nSamples;i++)
        {
            listWaveforms[Nwaveforms]->SetBinContent(i,buffer[0]);
            //totalOut << buffer[0] << endl;
            evtfile.read((char*)buffer,BufferBytes);
            //waveform[i]=buffer[0];
        }


        Nwaveforms++;
    }

    else
    {
        cout << "ERROR: unknown value for event type" << endl;
        return;
    }
}

int main()
{


    
    f = new TFile("eventTest.root","RECREATE");
    tree = new TTree("tree","");
    tree->SetAutoSave(0);

    outDPP = new TH1S("outDPP","outDPP",1024,0,60000);
    //TH1S* listWav[1000];

    tree->Branch("event",&de.timetag,"time/I:sgQ:lgQ:baseline");


    evtfile.clear();
    evtfile.open(name.c_str(),ios::binary);      

    if (evtfile.bad()) cout << "bad" << endl;
    if (evtfile.is_open()) cout << "open" << endl;
    if (evtfile.good()) cout << "good" << endl;
    if (evtfile.fail()) cout << "fail" << endl;

    evtfile.read((char*)buffer,BufferBytes);

    point = buffer;

    // start reading the evtfile for events
    while(!evtfile.eof())
    {
        int channelNum = readHeader(evtfile);
        cout << channelNum << endl;
        headerPrint(totalOut);
        timeOut << timetag << endl;

        // funnel event data differently based on channel

        switch (channelNum)
        {
            case 0: 
                // target-changer data
                headerPrint(targetChangerOut);
                unpack(evtfile, targetChangerOut);
                break;
            case 1:
                break;
            case 2:
                // monitor data
                headerPrint(monitorOut);
                unpack(evtfile, monitorOut);
                break;
            case 3:
                break;
            case 4:
                // detector data
                headerPrint(detectorOut);
                unpack(evtfile, detectorOut);
                break;
            case 5:
                break;
            case 6:
                // detector data
                headerPrint(detectorOut);
                unpack(evtfile, detectorOut);
                break;
            case 7:
                break;
            default:
                cout << "ERROR: unknown value for channel type" << endl;
                return 0;
                break;
        }

        /*
           for (int i=0;i<10;i++)
           {
           cout << *point << endl;
           point++;
           }
           */

        Ne++;
    }
    cout << "Finished processing buffer" << endl;
    f->Write();
    tree->Write();
    cout << Ne++ << endl;
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


