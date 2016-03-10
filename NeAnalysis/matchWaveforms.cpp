#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TROOT.h"
#include "TApplication.h"

using namespace std;

const string analysispath =  "/media/Drive3/";


/* event variables */

// variables for holding tree event data to be histogrammed
unsigned int evtNo, macroNo;
unsigned int sgQ, lgQ;
double completeTime;

vector<int> *waveform; // for holding one event's waveform data


// Re-link to an already-existing tree's data so we can read the tree
void setBranches(TTree* tree)
{
    tree->SetBranchAddress("macroNo",&macroNo);
    tree->SetBranchAddress("macroTime",&macroTime);
    //tree->SetBranchAddress("microNo",&microNo);
    tree->SetBranchAddress("evtNo",&evtNo);
    tree->SetBranchAddress("completeTime",&completeTime);
    //tree->SetBranchAddress("microTime",&microTime);
    tree->SetBranchAddress("sgQ",&sgQ);
    tree->SetBranchAddress("lgQ",&lgQ);
    tree->SetBranchAddress("waveform",&waveform);
}

void matchWaveforms()
{
    // we want to pull events from the overlap window where both ch4 and ch6
    // see the same events
    setBranches(ch6Tree);
    int totalEntries = ch6Tree->GetEntries();
    cout << "Searching through ch6Tree for events in the overlap region" << endl;

    // loop through the channel-specific tree and populate histos
    for(int j=0; j<totalEntries; j++)
    {
        ch6Tree->GetEntry(j);

        double timeDiff = completeTime-macroTime+TIME_OFFSET;

        double microTime = fmod(timeDiff,MICRO_PERIOD);

        if(microTime
        // if waveform data for this event exist, we want to populate
            // a histogram to display it

            // only plot 1 out of 10000 waveforms to save space and processing
            // time
            if(waveform->size() > 0 && j%10000 == 0)
            {
                waveformsDir = (TDirectory*)gDirectory->Get("waveformsDir");
                waveformsDir->cd();

                stringstream temp;
                temp << "macroNo " << macroNo << ", evtNo " << evtNo;
                waveformH = new TH1I(temp.str().c_str(),temp.str().c_str(),waveform->size(),0,waveform->size()*2);

                // loop through waveform data and fill histo
                for(int k=0; k<waveform->size(); k++)
                {
                    waveformH->SetBinContent(k,waveform->at(k));
                }
                gDirectory->cd("..");
            }
        }
    }

    // fill TOF, cross-section, etc. histos for channels 2, 4, 6
    for(int i=1; i<4; i++)
    {
        fillAdvancedHistos(i);
    }

    // fill basic histograms for waveform mode in each channel

    // first loop through all channel-specific waveform-mode trees
    for(int i=0; i<orchardW.size(); i++)
    {
        // create a channel-specific directory for each tree
        gDirectory->cd("/");
        dirs[i] = dirs[i] + "WaveformMode";
        gDirectory->mkdir(dirs[i].c_str(),dirs[i].c_str());

        gDirectory->GetDirectory(dirs[i].c_str())->cd();

        // instantiate histograms inside the channel-specific directory
        TH1I* macroNoH = new TH1I("macroNoH","macroNo",100000,0,100000);
        macroNoH->GetXaxis()->SetTitle("macropulse number of each event");

        TH1I* evtNoH = new TH1I("evtNoH","evtNo",2500,0,2500);
        evtNoH->GetXaxis()->SetTitle("event number of each event");

        TH1I* completeTimeH = new TH1I("completeTimeH","completeTime",6000,0,6000000000);
        completeTimeH->GetXaxis()->SetTitle("complete time for each event");

        // create subdirectory for holding waveform-mode waveform data
        gDirectory->mkdir("waveformsDir","raw DPP waveforms");
        waveformsDir = (TDirectory*)gDirectory->Get("waveformsDir");
        waveformsDir->cd();

        // fill waveform mode histograms
        setBranchesW(orchardW[i]);
        int totalEntries = orchardW[i]->GetEntries();
        cout << "Populating " << dirs[i] << " histograms..." << endl;

        // create a holder for the time-zero of each macropulse waveform
        // so we can plot all 11 waveform chunks on the same histogram
        // initialize to -650000 to ensure that the first event in the waveform
        // trees is always considered the start of a macropulse waveform
        double waveformStart = -650000;

        // we need label the number of waveform-mode macropulses to make
        // uniquely named histograms
        int waveformNo = 0;

        for(int j=0; j<totalEntries; j++)
        {
            //cout << "looping thru waveform mode events, #" << j << endl;
            //fflush(stdout);

            orchardW[i]->GetEntry(j);

            macroNoH->Fill(macroNo);
            evtNoH->Fill(evtNo);
            completeTimeH->Fill(completeTime);

            if(completeTime >= waveformStart+650000)
            {
                // new macropulse in waveform mode because no other waveform
                // data came before it in the last macropulse period
                // create a new plot to stitch 11 waveforms together
                stringstream temp;
                temp << "full waveform " << waveformNo;
                waveformH = new TH1I(temp.str().c_str(),temp.str().c_str(),350000,0,700000);

                // set the start of the macropulse to the first event timer
                waveformStart = completeTime;
                waveformNo++;
            }

            for(int k=0; k<waveform->size(); k++)
            {
                waveformH->SetBinContent(k+(completeTime-waveformStart)/2,waveform->at(k));
            }
        }
    }
}

int main(int argc, char* argv[])
{
    // read in the raw file name
    string runDir = argv[1];
    string runNo = argv[2];

    // Open the raw tree from the initial sort. If it doesn't exist, exit.
    stringstream treeName;
    stringstream fileInName;
    stringstream fileOutName;
    stringstream scavengerEventsName;
    stringstream summedDetEventsName;

    treeName << runDir << "-" << runNo; 

    fileInName << analysispath <<"analysis/" << runDir << "/" << treeName.str() << "_sorted.root";
    fileOutName << analysispath <<"analysis/" << runDir << "/" << treeName.str() << "_ch4ch6_matchedWaveforms.root";

    TFile* file = new TFile(fileInName.str().c_str(),"READ");

    if(file->Get("ch4Tree") && file->Get("ch6Tree"))
    {
        cout << "Located ch4/6 trees in " << fileInName << "." << endl;

    TTree* ch4Tree = (TTree*)file->Get("ch4Tree");
    TTree* ch6Tree = (TTree*)file->Get("ch6Tree");
    }

    else
    {
        cout << "Failed to find ch4/6 trees in " << fileInName << ". Exiting." << endl;
        exit(1);
    }
   
    // open output file to contain histos
    TFile* fileOut = new TFile(fileOutName.str().c_str(),"RECREATE");

    // pull and histogram the first 50 events from ch4 and ch6 that are in both
    // channels (i.e., the same event is recorded in both channel 4 and channel
    // 6 because it's in the small window where both channels overlap)
    matchWaveforms();

    fileOut->Write();

    fileOut->Close();
}
