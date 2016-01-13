// Takes a ROOT tree with raw event data and processes it into a new tree useful
// for analysis

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TEntryList.h"

using namespace std;

// Experimental constants


// Time delay of target changer coarse time after real macropulse start time
const double TIME_OFFSET = 836; // in ns

const string analysispath =  "/media/Drive3/";


// variables for holding raw tree event data
unsigned int chNo, evtType, extTime, fineTime, sgQ, lgQ;
double timetag;

// additional variables for holding new tree event data
unsigned int evtNo, macroNo, targetPos;
double completeTime, macroTime;

vector<int> waveform; // for holding one event's waveform data

// event structure for resorted trees (different from raw tree event structure)
struct event
{
  unsigned int macroNo; // label each event by the macropulse it's in
  unsigned int evtNo; // uniquely label each event in a macropulse

  double macroTime; // =extTime+timetag for target changer
  double completeTime; // =extTime+timetag+fineTime for detector events

  unsigned int targetPos; // target changer position;

  unsigned int sgQ, lgQ; // event charge gates

  vector<int> waveform;
} ev;

TTree* cTree; // holds raw data with beam-off periods removed

// channel-specific trees for DPP and waveform mode
// data from cTree will be reprocessed and divided into these trees for later
// analysis/histogramming
TTree* ch0Tree;
TTree* ch2Tree;
TTree* ch4Tree;
TTree* ch6Tree;
TTree* ch0TreeW;
TTree* ch2TreeW;
TTree* ch4TreeW;
TTree* ch6TreeW;


// Target Changer lgQ gates for setting integral target positions 1-6
const int tarGate[12] = {5000,7500,11000,13500,17000,20000,23000,27000,29000,33000,35000,39000};

void fillTree(TTree* tree)
{
    ev.macroNo = macroNo;
    ev.evtNo = evtNo;
    ev.macroTime = macroTime;
    ev.completeTime = completeTime;
    ev.targetPos = targetPos;
    ev.sgQ = sgQ;
    ev.lgQ = lgQ;
    ev.waveform = waveform; 

    tree->Fill();
}

void cleanTree(TTree* tree)
{
    // This function excises periods of beam off from the raw tree
    // The approach is:
    //  - to loop through all entries in the tree in the order that
    //    they were read out of the digitizer
    //  - when a target changer event is found that doesn't line up
    //    with the macropulse structure (i.e., timestamps 8.3 ms or
    //    16.6 ms after each other), this indicates that beam went
    //    off in the previous macropulse.
    //  - Move backward in the tree until the previous chunk of
    //    target changer events is found, and save the event indices
    //    of the beginning and end of this dead time
    //  - Using this info, fill a vector with the periods when the
    //    beam is on, and recreate the raw tree using these periods
    //    as gates for when data should be allowed. The resulting
    //    tree, "cTree", is used for all analysis.

    // link raw tree branches to variables-to-be-read
    tree->SetBranchAddress("chNo",&chNo);
    tree->SetBranchAddress("sgQ",&sgQ);
    tree->SetBranchAddress("lgQ",&lgQ);
    tree->SetBranchAddress("evtType",&evtType);
    tree->SetBranchAddress("timetag",&timetag);
    tree->SetBranchAddress("extTime",&extTime);
    tree->SetBranchAddress("fineTime",&fineTime);
    tree->SetBranchAddress("waveform",&waveform);

    // prepare variables for evaluating whether a macropulse has been skipped
    // (thus indicating that beam went off)
    double timeDiff = 0; // difference between times of adjacent macropulses
    double prevMacroTime = 0; // placeholder for previous macropulse time
    double currentMacroTime = 0; // holds complete timestamp of current macro
    int init = 0; // for holding the first event index of a 'beam on' period
    int prevEvtType = 2; // keeps track of when acquisition switches modes;
                         // needed to accept the first macropulse of each
                         // new acquisition pariod
    vector<pair<int,int>> beamOn; // holds event indices for when beam is on
                         // first index is 'start of beam on', second index
                         // is 'end of beam on'

    // loop through raw tree looking for target changer events
    int totalEntries = tree->GetEntries();

    for(int i = 0; i<totalEntries; i++)
    {
        tree->GetEntry(i);

        if (chNo == 0)
        {
            // found a target changer event - now test for beam on/off
            currentMacroTime = ((double)extTime*pow(2,32)+timetag);
            timeDiff = currentMacroTime-prevMacroTime;

            if ((timeDiff > 8250000 && timeDiff < 8350000) || (timeDiff > 16550000 && timeDiff < 16650000) || prevEvtType == 2)
            {
                // the target changer event came within the expected window of
                // 8.3 or 16.6 ms, or was the start of a new acquisition period
                // ...thus beam was on, so accept this event and continue
            }

            else
            {
                // target changer time was OUTSIDE acceptable bounds
                // thus beam was off for some period, so figure out that period

                // we'll need to go backwards to find the previous set of
                // target changer events when beam was still good
                // shift the event index backwards by 10 so we can escape the current
                // chunk of target changer events where beam off was detected
                for (int j = i-10; j>0; j--)
                {
                    tree->GetEntry(j);
                    if(chNo=0)
                    {
                        // found the previous chunk of target changer events;
                        // record their index in beamOn
                        beamOn.push_back(make_pair(init,j));

                        // shift init to the next event in preparation for next
                        // beam on period
                        init = i+1;

                        // ... and keep looping forward through raw
                        // tree where we left off (that is, at index i)
                        break;
                    }
                }
            }

            // reset placeholder variables for next event
            prevEvtType = evtType;
            prevMacroTime = currentMacroTime;

        }
    }

    // after loop finishes, add the last leg of the run to the beamOn list
    beamOn.push_back(make_pair(init,totalEntries));

    // create a new empty tree to hold the cleaned data from the raw tree
    cTree = tree->CloneTree(0);

    // Add only events when beam was on to the cleaned tree
    for(int i = 0; i<beamOn.size(); i++)
    {
        int resume = beamOn[i].first;
        int omit = beamOn[i].second;

        for(int j = beamOn[i].first; j<beamOn[i].second; j++)
        {
            tree->GetEntry(j);
            cTree->Fill();
        }
    }
}

void populateTrees()
{
    // This method creates two trees for each channel, one for DPP mode and one
    // for waveform mode. Each event to be added to these trees will be
    // assigned to a macropulse for easier analysis later.

    // create subsets of the cleaned tree based on channel number
    cTree->Draw(">>targetChEvents","chNo==0","entrylist");
    cTree->Draw(">>monitorEvents","chNo==2","entrylist");
    cTree->Draw(">>detSEvents","chNo==4","entrylist");
    cTree->Draw(">>scavengerEvents","chNo==6","entrylist");

    TEntryList* targetChEvents = (TEntryList*)gDirectory->Get("targetChEvents");
    TEntryList* monitorEvents = (TEntryList*)gDirectory->Get("monitorEvents");
    TEntryList* detSEvents = (TEntryList*)gDirectory->Get("detSEvents");
    TEntryList* scavengerEvents = (TEntryList*)gDirectory->Get("scavengerEvents");

    // add these event subsets for a list for easier looping
    vector<TEntryList*> channelList;
    channelList.push_back(targetChEvents);
    channelList.push_back(monitorEvents);
    channelList.push_back(detSEvents);
    channelList.push_back(scavengerEvents);

    for(int j = 0; j<channelList.size(); j++)
    {
        // looping through tree once per channel number
        cout << "Populating trees" << endl;

        cTree->SetEntryList(channelList[j]);
        int channelEntries = cTree->GetEntries();

        int prevEvtType = 0; // keeps track of when acquisition switches modes
        // so that macropulse counting is done correctly

        int prevTime = 0;         // keep track of previous macro time 

        switch (j)
        {
            case 0:
                // Target changer (ch0)
                // Populate a tree of only target changer events in preparation
                // for assigning times and target changer positions to the
                // detector channels (ch2, ch4, ch6)

                for (int i=0; i<channelEntries; i++)
                {
                    cTree->GetEntry(i);

                    completeTime = (double)extTime*pow(2,32)+timetag;
                    macroTime = completeTime;

                    if (evtType==1)
                    {
                        // new macropulse in DPP mode

                        // assign an integral target position based on lgQ value
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
                            targetPos = 0;
                        }

                        if (prevEvtType==2)
                        {
                            // the previous event was a waveform event, so this
                            // is the first DPP mode event

                            // clear the evtNo counter
                            evtNo = 0;
                        }

                        fillTree(ch0Tree);

                        macroNo++; // increment macropulse counter
                    }

                    else if (evtType==2)
                    {
                        // waveform mode
                        // don't increment the macropulse number, but fill
                        // the tree with timestamp, evtNo, and waveform data
                        targetPos = 0;

                        if (prevEvtType==1)
                        {
                            // the previous event was a DPP event, so this is
                            // the first waveform mode event in the wavelet

                            // clear the evtNo counter
                            evtNo = 0;
                            prevTime = timetag+pow(2,32)*extTime;
                        }

                        else if (timetag+pow(2,32)*extTime>prevTime+7000000)
                        {
                            evtNo = 0; // clear the evtNo counter
                        }

                        fillTree(ch0TreeW);
                    }

                    evtNo++;

                    prevEvtType = evtType;

                    if (macroNo%100==0)
                    {
                        cout << "Macro number " << macroNo << ", fullTime " << timetag+pow(2,32)*extTime << "\r";
                        fflush(stdout);
                    }
                }

                cout << endl;

                break;

            case 1:
                // Monitor (ch2)
                // Populate a tree of only monitor events

                // link tree branches to variables-to-read
                ch0Tree->SetBranchAddress("macroNo",&macroNo);
                ch0Tree->SetBranchAddress("completeTime",&macroTime);
                ch0Tree->SetBranchAddress("targetPos",&targetPos);

                ch0Tree->GetEntry(0);
                // we're pointed at the first macropulse in preparation for
                // scanning through all monitor events

                for (int i=0; i<channelEntries; i++)
                {
                    cTree->GetEntry(i);

                    completeTime = (double)extTime*pow(2,32)+timetag;

                    if (evtType==1)
                    {
                        while (completeTime-macroTime+TIME_OFFSET > 8000000)
                        {
                            // previous macropulse has elapsed -
                            // move to the next target changer event
                            ch0Tree->GetEntry(macroNo+1);
                        }

                        fillTree(ch2Tree);
                    }

                    else if (evtType==2)
                    {
                        if (prevEvtType==1)
                        {
                            ch0Tree->GetEntry(macroNo+1); // shift to the next
                            // macropulse in preparation for the next DPP mode
                        }

                        fillTree(ch2TreeW);
                    }

                    evtNo++; // increment event counter for channel 2

                    prevEvtType = evtType;
                }

                break;

            case 2: // summed detector (ch4)
            case 3: // scavenger (ch6)

                // link tree branches to variables-to-read
                ch0Tree->SetBranchAddress("macroNo",&macroNo);
                ch0Tree->SetBranchAddress("completeTime",&macroTime);
                ch0Tree->SetBranchAddress("targetPos",&targetPos);

                ch0Tree->GetEntry(0);
                // we're pointed at the first macropulse in preparation for
                // scanning through all detector events

                for (int i=0; i<channelEntries; i++)
                {
                    cTree->GetEntry(i);

                    completeTime = (double)extTime*pow(2,32)+timetag+(double)fineTime*(2./1024.);

                    if (evtType==1)
                    {
                        while (completeTime-macroTime+TIME_OFFSET > 8000000)
                        {
                            // previous macropulse has elapsed -
                            // move to the next target changer event
                            ch0Tree->GetEntry(macroNo+1);
                        }

                        if (chNo==4)
                        {
                            fillTree(ch4Tree);
                        }

                        else
                        {
                            fillTree(ch6Tree);
                        }
                    }

                    else if (evtType==2)
                    {
                        if (prevEvtType==1)
                        {
                            ch0Tree->GetEntry(macroNo+1); // shift to the next
                            // macropulse in preparation for the next DPP mode
                        }

                        if (chNo==4)
                        {
                            fillTree(ch4TreeW);
                        }

                        else
                        {
                            fillTree(ch6TreeW);
                        }

                        evtNo++; // increment event counter for channel 2

                        prevEvtType = evtType;
                    }

                    break;

                }

            default:
                break;
        }
    }
}

void branch(TTree* tree)
{
    tree->Branch("macroNo",&ev.macroNo,"macroNo/i");
    tree->Branch("evtNo",&ev.evtNo,"evtNo/i");
    tree->Branch("macroTime",&ev.macroTime,"macroTime/i");
    tree->Branch("completeTime",&ev.completeTime,"completeTime/d");
    tree->Branch("targetPos",&ev.targetPos,"targetPos/i");
    tree->Branch("sgQ",&ev.sgQ,"sgQ/i");
    tree->Branch("lgQ",&ev.lgQ,"lgQ/i");
    tree->Branch("waveform",&ev.waveform);
}


int main(int argc, char* argv[])
{

    // Open the raw tree from an initial sort. If it doesn't exist, exit.
    TFile *file;
    TTree *tree;

    string runDir, runNo;
    stringstream treeName;
    stringstream fileName;

    treeName << runDir << "-" << runNo; 
    fileName << analysispath <<"analysis/" << runDir << "/" << treeName.str() << ".root";

    file = new TFile(fileName.str().c_str(),"UPDATE");

    if(file->Get("tree"))
    {
        // Found the raw tree; start sorting the raw trees into new trees for
        // each channel.
        cout << "Found previous raw tree; creating processed trees for " << fileName << endl;

        tree = (TTree*)file->Get("tree");

        ch0Tree = new TTree("ch0Tree","");
        branch(ch0Tree);

        ch2Tree = new TTree("ch2Tree","");
        branch(ch2Tree);

        ch4Tree = new TTree("ch4Tree","");
        branch(ch4Tree);

        ch6Tree = new TTree("ch6Tree","");
        branch(ch6Tree);

        ch0TreeW = new TTree("ch0TreeW","");
        branch(ch0TreeW);

        ch2TreeW = new TTree("ch2TreeW","");
        branch(ch2TreeW);

        ch4TreeW = new TTree("ch4TreeW","");
        branch(ch4TreeW);

        ch6TreeW = new TTree("ch6TreeW","");
        branch(ch6TreeW);

    }

    else
    {
        cout << "Failed to find raw tree - exiting " << fileName << endl;
        exit(1);
    }

    cleanTree(tree);

    populateTrees();

    file->Write();

    file->Close();
}
