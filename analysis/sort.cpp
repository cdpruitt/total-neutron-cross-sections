// C++ file to read CAEN event files 

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TH1S.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

int main()
{

    TFile* f = new TFile("eventTest.root","RECREATE");
    TTree* tree = new TTree("tree","");
    tree->SetAutoSave(0);

    TH1S* outDPP = new TH1S("outDPP","outDPP",1024,0,60000);
    //TH1S* listWav[1000];
    vector<TH1S*> listWaveforms;
    vector<TH1S*> listWavelets; 

    struct DPPevent {
        unsigned int timetag, sgQ, lgQ, baseline;
    };

    DPPevent de;

    tree->Branch("event",&de.timetag,"time/I:sgQ:lgQ:baseline");

    int const BufferWords = 1;
    int const BufferBytes = BufferWords*2;
    unsigned short buffer[BufferWords];
    unsigned short *point;

    ostringstream histName;
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

    unsigned int Ne = 0; // counter for the number of events in the buffer
    unsigned int Nwavelets = 0; // counter for the number of MIXED mode wavelets in the buffer
    unsigned int Nwaveforms = 0; // counter for the number of waveform events in the buffer

    while(!evtfile.eof())
    {
        // start reading header

        // size is the number of bytes in the event (self-inclusive)
        unsigned short size1 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        unsigned short size2 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        unsigned long size = (size2 << 16) | size1;
        cout << "size = " << size << endl;

        // evtype is either 1 (DPP data) or 2 (waveform data)
        unsigned short evtype1 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        unsigned short evtype2 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        unsigned long evtype = (evtype2 << 16) | evtype1;
        cout << "event type = " << evtype << endl;

        // channel ranges from 0-7; each channel is considered an independent event generator
        unsigned short channel1 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        unsigned short channel2 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        unsigned long channel = (channel2 << 16) | channel1;
        cout << "channel type = " << channel << endl;

        // timetag is the coarse trigger time for this event, in 2ns increments
        unsigned short timetag1 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        unsigned short timetag2 = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);
        unsigned long timetag = (timetag2<< 16) | timetag1;
        cout << "timetag = " << timetag << endl;

        if(evtype==1)
        {
            cout << "DPP event data" << endl;

            // zcTime is the interpolated zero-crossing time, in picoseconds
            unsigned short zcTime = buffer[0];
            evtfile.read((char*)buffer,BufferBytes);
            cout << "zero crossing Time = " << zcTime << endl;

            // sgQ is the short gate integrated charge, in digitizer units 
            unsigned short sgQ = buffer[0];
            evtfile.read((char*)buffer,BufferBytes);
            cout << "short gate charge = " << sgQ << endl;

            // lgQ is the long gate integrated charge, in digitizer units 
            unsigned short lgQ = buffer[0];
            evtfile.read((char*)buffer,BufferBytes);
            outDPP->Fill(lgQ);
            cout << "long gate charge = " << lgQ << endl;

            // baseline is the baseline level, frozen at the trigger time
            unsigned short baseline = buffer[0];
            evtfile.read((char*)buffer,BufferBytes);
            cout << "baseline level = " << baseline << endl;

            // puRej is a pile-up rejection flag (1 if pile-up detected? 0 if not?)
            unsigned short puRej = buffer[0];
            evtfile.read((char*)buffer,BufferBytes);
            cout << "pile-up detected = " << puRej << endl;

            // nSamp is the number of waveform samples that follow (in LIST mode, this is 0)
            unsigned short nSamp1 = buffer[0];
            evtfile.read((char*)buffer,BufferBytes);
            unsigned short nSamp2 = buffer[0];
            evtfile.read((char*)buffer,BufferBytes);
            const unsigned int nSamp = (nSamp2 << 16) | nSamp1;
            cout << "number of waveform samples = " << nSamp << endl;

            if(nSamp > 0)
                // We must be in MIXED mode; time to output wavelet.
            {
                cout << "MIXED event wavelet data" << endl;

                histName.str("");
                histName << "outWavelet" << Nwavelets;
                string tempHist = histName.str();

                listWavelets.push_back(new TH1S(tempHist.c_str(),tempHist.c_str(),nSamp,0,nSamp*2));

                for(int i=0;i<nSamp;i++)
                {
                    listWavelets[Nwavelets]->SetBinContent(i,buffer[0]);
                    cout << buffer[0] << endl;
                    evtfile.read((char*)buffer,BufferBytes);
                    //waveform[i]=buffer[0];
                }
                Nwavelets++;
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

            de.timetag = timetag;
            de.sgQ = sgQ;
            de.lgQ = lgQ;
            de.baseline = baseline;

            tree->Fill();

        }

        else if(evtype==2)
        {
            cout << "Waveform event data" << endl;

            unsigned short nSamples1 = buffer[0];
            evtfile.read((char*)buffer,BufferBytes);
            unsigned short nSamples2 = buffer[0];
            evtfile.read((char*)buffer,BufferBytes);
            unsigned long nSamples = (nSamples2 << 16) | nSamples1;
            cout << "nSamples = " << nSamples << endl;

            //unsigned short waveform[nSamples];

            histName.str("");
            histName << "outWaveform" << Nwaveforms;
            string tempHist = histName.str();

            listWaveforms.push_back(new TH1S(tempHist.c_str(),tempHist.c_str(),nSamples,0,nSamples*2));

            for(int i=0;i<nSamples;i++)
            {
                listWaveforms[Nwaveforms]->SetBinContent(i,buffer[0]);
                cout << buffer[0] << endl;
                evtfile.read((char*)buffer,BufferBytes);
                //waveform[i]=buffer[0];
            }


            Nwaveforms++;
        }

        else
        {
            cout << "ERROR: unknown value for event type" << endl;
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
        cout << endl;
    }
    cout << "Finished processing buffer" << endl;
    f->Write();
    tree->Write();
    cout << Ne++ << endl;
}
