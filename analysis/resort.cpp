// Takes a ROOT tree with raw event data and processes it into channel-specific
// trees useful for analysis

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TEntryList.h"

using namespace std;

// file path to event data
const string analysispath =  "/media/Drive3/";



// EXPERIMENTAL PARAMETERS

// Time delay of target changer coarse time after real macropulse start time
const double TIME_OFFSET = 836; // in ns

// Period of micropulses
const double MICRO_PERIOD = 1788.820; // in ns

// Period of macropulses
const double MACRO_PERIOD = 8330000; // in ns

// Target Changer lgQ gates for setting integral target positions 1-6
const int tarGate[12] = {5000,7500,11000,13500,17000,20000,23000,27000,29000,33000,35000,39000};



// TREE VARIABLES

// for holding event data from the raw tree
unsigned int chNo, evtType, extTime, fineTime, sgQ, lgQ;
double timetag;

// for holding event data for the new channel-specific trees
unsigned int evtNo, macroNo, microNo, targetPos;
double completeTime, macroTime, microTime;
vector<int> waveform;

// channel-specific trees' event structure (different from raw tree's)
struct event
{
  unsigned int macroNo; // label each event by macropulse
  unsigned int microNo; // label each event by micropulse
  unsigned int evtNo; // uniquely label each event in a macropulse

  double macroTime; // the event's time-zero reference (the macropulse start)
  double completeTime; // the event's 48-bit timestamp
  double microTime; // the event's timestamp relative to the micropulse it's in

  unsigned int targetPos; // target position

  unsigned int sgQ, lgQ; // the event's short and long integrated charge gates

  vector<int> waveform; // waveform data for this event
} ev;

// for transferring waveform data from clean tree to subtrees
vector<int> *dummyWaveform;

// to assign the correct macropulses to each event, we need to track timestamp
// resets (DPP/waveform mode changes). Use these 'holder' variables to compare
// consecutive detector and target changer timestamps.
double prevMacroTime;
double prevCompleteTime;
int prevEvtType;

// detector index (2 is monitor, 4 is summed detector, 6 is scavenger)
int detIndex;

// TREES

// raw tree
// tree points at the already-populated TTree 'tree' that was filled in raw.cpp
// it contains unprocessed events straight from the digitizer's binary output
// by looping through this tree, we'll segregate each event by channel to
// perform later analysis (assign TOFs, calculate cross-sections, etc)
TTree *tree;

// channel-specific trees
// Each channel has a tree for DPP and waveform mode (except for channel 6,
// which collected no waveform mode data)
TTree* ch0Tree;
TTree* ch2Tree;
TTree* ch4Tree;
TTree* ch6Tree;
TTree* ch0TreeW;
TTree* ch2TreeW;
TTree* ch4TreeW;


/************ METHODS *************/

// connect a channel-specific tree to correct variables so we can populate it
void branch(TTree* tree)
{
    tree->Branch("macroNo",&ev.macroNo,"macroNo/i");
    tree->Branch("microNo",&ev.microNo,"microNo/i");
    tree->Branch("evtNo",&ev.evtNo,"evtNo/i");
    tree->Branch("macroTime",&ev.macroTime,"macroTime/d");
    tree->Branch("completeTime",&ev.completeTime,"completeTime/d");
    tree->Branch("microTime",&ev.microTime,"microTime/d");
    tree->Branch("targetPos",&ev.targetPos,"targetPos/i");
    tree->Branch("sgQ",&ev.sgQ,"sgQ/i");
    tree->Branch("lgQ",&ev.lgQ,"lgQ/i");
    tree->Branch("waveform",&ev.waveform);
}

// connect a channel-specific tree to correct variables so we can populate it
void branchW(TTree* tree)
{
    tree->Branch("macroNo",&ev.macroNo,"macroNo/i");
    tree->Branch("evtNo",&ev.evtNo,"evtNo/i");
    tree->Branch("completeTime",&ev.completeTime,"completeTime/d");
    tree->Branch("waveform",&ev.waveform);
}

// add a new event to a channel-specific tree.
void fillTree(TTree* tree)
{
    ev.macroNo = macroNo;
    ev.microNo = microNo;
    ev.evtNo = evtNo;
    ev.macroTime = macroTime;
    ev.completeTime = completeTime;
    ev.microTime = microTime;
    ev.targetPos = targetPos;
    ev.sgQ = sgQ;
    ev.lgQ = lgQ;

    /*if ((int)completeTime%10000 != 0 && evtType==1)
    {
        waveform.clear(); 
    }*/

    ev.waveform = waveform; 

    tree->Fill();
}

// pull an event's waveform data in preparation for storing it in the
// channel-specific trees
void getWaveform(vector<int>* dummyWaveform)
{
    // empty the stale waveform data
    waveform.clear();

    // fill waveform with wavelet data from an event
    // (the raw tree's branch points to dummyWaveform)
    for(int k = 0; k<dummyWaveform->size(); k++)
    {
        waveform.push_back(dummyWaveform->at(k));
    }
}

// add a waveform-mode event to a channel specific tree
void processDPPEvent(TTree* tree)
{
    // pointing at the correct macropulse - fill the
    // channel-specific tree with data from this channel

    // calculate the time elapsed since the start of the
    // current macropulse and the time elapsed since the start
    // of the micropulse
    double timeDiff = completeTime-macroTime-TIME_OFFSET;

    microTime = fmod(timeDiff,MICRO_PERIOD);
    microNo = floor(timeDiff/MICRO_PERIOD);

    getWaveform(dummyWaveform);

    fillTree(tree);
}

// add a waveform-mode event to a channel specific tree
void processWaveformEvent(TTree* tree)
{
    // waveform mode

    // update waveform data
    getWaveform(dummyWaveform);

    if (prevEvtType==1)
    {
        evtNo=0;
    }

    fillTree(tree);
}



// using the lgQ from the target changer, assign each macropulse a target
// position
int assignTargetPos(int lgQ)
{
    // assign an integral target position based on lgQ value
    // passed in by a target changer event
    for(int k = 0; k < 6; k ++)
    {
        if (lgQ>tarGate[2*k] && lgQ<tarGate[2*k+1])
        {
            // lgQ fits within one of the lgQ gates we defined, so
            // assign that target position
            return k;
        }

        else
        {
            // lgQ doesn't fit within one of the lgQ gates, so set
            // the target position to 0. These will be thrown out
            // during plotting of cross-sections (no target was
            // in-beam during this event).
            return 0;
        }
    }
}

// fill the ch0Tree with target changer events so that we have macropulse
// start times and target positions available for later processing
// 
// this is the first channel-specific tree we fill
void processTargetChanger()
{
    // Populate a tree of only target changer events in preparation
    // for assigning times and target changer positions to the
    // detector channels (ch2, ch4, ch6)

    // set macroNo to 0 for the start of the run
    macroNo = 0;

    // loop through raw tree to find target changer events
    int totalEntries = tree->GetEntries();

    for (int i=0; i<totalEntries; i++)
    {
        tree->GetEntry(i);

        if (chNo == 0)
        {
            // found a target changer event
            // calculate the time of the start of the macropulse using the
            // target changer event
            macroTime = (double)extTime*pow(2,32)+timetag;
            completeTime = macroTime;

            // figure out which target is in the beam during this event
            targetPos = assignTargetPos(lgQ);

            // extract the waveform data for the event
            getWaveform(dummyWaveform);

            // check to see whether the event is DPP or waveform mode
            if (evtType==1)
            {
                // DPP mode

                // First, check to see if this is the start of a new DPP period
                // (i.e., a switch between DPP and waveform modes)
                if(prevEvtType==2)
                {
                    // the last event from this channel was in waveform mode
                    // so we must be in a new DPP mode period.
                    // Set evtNo = 1 to indicate this.
                    // When we sort through the other channel specific
                    // trees we'll know where the mode switches are. 
                    evtNo = 1;
                }

                else
                {
                    // NOT the start of a new DPP period (99.8% of the time)
                    evtNo = 0;
                }

                // all the variables are updated, so fill the target changer
                // tree with the event
                fillTree(ch0Tree);

                // increment macropulse counter to prepare for next macropulse
                macroNo++;
            }

            else if (evtType==2)
            {
                // waveform mode

                // Use the waveform-mode data to fill a waveform-only
                // tree (ch0TreeW).

                // To plot a macropulse-long waveform, we'll  need to stitch
                // together several consecutive wavelets.
                // The time-zero for this stitched-together macropulse waveform
                // should be from the first event in the train of wavelets.
                double zeroTime;

                // First we must check if this is either:
                //  - the first wavelet in a train of wavelets, or
                //  - at least 700 us have elapsed since the last macropulse, so
                //    we're in a new macropulse
                if (prevEvtType==1 || macroTime>prevMacroTime+7000000)
                {
                    // either just started waveform mode or this is a new
                    // macropulse; reset the event number
                    evtNo = 0;
                }

                else
                {
                    // this is NOT the first wavelet in a wavelet train;
                    // increment the evtNo and proceed to fill the tree
                    evtNo++;
                }

                fillTree(ch0TreeW);
            }

            // before we pull the next event, save this event's type so we can
            // check to see if the new event changes DPP/waveform mode
            prevEvtType = evtType;
        }
        // move to next event in the loop
    }
    // loop finished
    // all target changer events sorted into ch0Tree and ch0TreeW
    // time to move to detector channels
}

// loop over the raw tree containing all events, and assign detector events
// (as opposed to target changer events) to the correct channel-specific tree
void processDetEvents()
{
    // Populate events from detector index 2, 4, 6 (monitor, summed detector,
    // and scavenger, respectively) into channel-specific trees.

    for(int detIndex = 2; detIndex < 8; detIndex+=2)
    {
        // To prepare detector events for later analysis, we'll determine which
        // macropulse and micropulse each event is in.
        // The time reference for each detector event will be the macropulse start
        // time (macroTime).

        // each event is uniquely labeled by its channel-specific order in the
        // macropulse
        // reset event counter in preparation for looping through events
        evtNo = 0; 

        // point to the first macropulse in preparation for reading detector events
        ch0Tree->GetEntry(0);

        // reset the timestamp holders
        prevMacroTime = 0;
        prevCompleteTime = 0;

        TTree *chTree, *chTreeW;
        switch(detIndex)
        {
            case 2:
                chTree = ch2Tree;
                chTreeW = ch2TreeW;
                break;
            case 4:
                chTree = ch4Tree;
                chTreeW = ch4TreeW;
                break;
            case 6:
                chTree = ch6Tree;
                // no ch6 waveform tree
                break;
        }

        // all the preliminary steps ready:
        // start looping through detector events
        int totalEntries = tree->GetEntries();

        for (int index = 0; index<totalEntries; index++)
        {
            // pull next event's data from the raw tree
            tree->GetEntry(index);

            if (chNo == detIndex)
            {
                // we've located a detector event from the correct (jth) channel.
                // use its data along with the target changer data to assign this
                // event the correct time, macropulse no, etc.

                // first calculate this event's time
                completeTime = (double)extTime*pow(2,32)+timetag;

                // now check to see if the event is DPP or waveform mode
                if (evtType == 1)
                {
                    // DPP mode

                    // first, check to see if we've wrapped around
                    // to a new DPP mode period.
                    // when a DPP/waveform mode switch occurs, the timestamps are reset to 0 at
                    // the start of the new mode. To keep track of when a channel switches
                    // modes, we need to compare consecutive events' timestamps.
                    if (completeTime < prevCompleteTime)
                    {
                        // new DPP mode period, triggered by a
                        // channel-specific event wrapping back around
                        // to time zero

                        cout << "channel timestamp wraparound, completeTime = " << completeTime << ", prevCompleteTime = " << prevCompleteTime << endl;

                        // update the ch0Tree to point at the new
                        // macropulse
                        do
                        {
                            prevMacroTime = macroTime;
                            ch0Tree->GetEntry(macroNo+1);
                        }
                        while (prevMacroTime < macroTime);

                        cout << "start new DPP mode period: macroTime = " << macroTime << ", macroNo = " << macroNo << endl;
                        // ch0Tree is now pointing at the first
                        // macropulse of the next DPP mode period

                        // reset the event counter for this channel
                        evtNo = 0;
                    }

                    // check to see if enough time has elapsed such that
                    // we should update the current macropulse to a new
                    // one
                    else if (completeTime-macroTime-TIME_OFFSET > 8000000)
                    {
                        // macropulse expired
                        // time to examine the next macropulse

                        // hold the old macro time so we'll be able to
                        // compare it with the new macro time and make
                        // sure there's still normal beam structure
                        prevMacroTime = macroTime;

                        // pull the next macropulse 
                        ch0Tree->GetEntry(macroNo+1);
                        cout << "macroTime = " << macroTime << ", macroNo = " << macroNo << "\r";
                        fflush(stdout);

                        // first, we should check to make sure we're not
                        // about to run off the end of ch0Tree and get
                        // to the end of the run
                        if(macroNo+1>=ch0Tree->GetEntries())
                        {
                            // reached the end of the ch0Tree
                            // break to throw away all subsequent events in this
                            // channel because we've reached the end of the run
                            break;
                        }

                        // we're not at the end of the run, so we
                        // should examine the next macropulse to see
                        // what to do with it

                        // evtNo was updated when we pulled the next
                        // macropulse - let's check to see if we've
                        // reached the end of a DPP period
                        if (evtNo==1)
                        {
                            // new DPP period found because the
                            // macropulse has evtNo==1, so move the
                            // cleaned tree indices forward to the new
                            // DPP period

                            cout << "ch0Tree wraparound, macroTime = " << macroTime << ", prevMacroTime = " << prevMacroTime << endl;

                            do
                            {
                                do
                                {
                                    index++;
                                    tree->GetEntry(index);
                                }
                                while (chNo != detIndex || evtType != 1);

                                prevCompleteTime = completeTime;
                                completeTime = (double)extTime*pow(2,32)+timetag;
                            }
                            while (prevCompleteTime < completeTime);

                            cout << "Start new DPP period, completeTime = " << completeTime << ", prevCompleteTime = " << prevCompleteTime << endl;
                        }

                        // not a new DPP period, so lastly let's check
                        // to see if there are any beam anomalies (that
                        // is, the new macropulse is out-of-sync with
                        // the old one)
                        else if ((macroTime-prevMacroTime > 8300000 && macroTime-prevMacroTime < 8360000) || (macroTime-prevMacroTime > 16600000 && macroTime-prevMacroTime < 17200000) || (macroTime-prevMacroTime > 24900000 && macroTime-prevMacroTime < 25080000))
                        {
                            // no beam anomaly - keep filling the
                            // channel-specific tree and continue as
                            // normal
                        }

                        else
                        {
                            // there IS a beam anomaly - cycle through
                            // the cleaned and ch0Tree until we find
                            // another period of clean beam
                            while(!(macroTime-prevMacroTime > 8300000 && macroTime-prevMacroTime < 8360000) && !(macroTime-prevMacroTime > 16600000 && macroTime-prevMacroTime < 17200000) && !(macroTime-prevMacroTime > 24900000 && macroTime-prevMacroTime < 25080000))
                            {
                                // update the ch0Tree to point at the
                                // next macropulse so we can check to
                                // see if it's out-of-sync with the
                                // previous macropulse 
                                prevMacroTime = macroTime;
                                ch0Tree->GetEntry(macroNo+1);

                                cout << "Beam anomaly detected." << endl;
                                cout << "macroTime is " << macroTime << " macroNo " << macroNo << endl;
                            };

                            // OK - ch0Tree is now in a period of
                            // good beam. We need to move the cleaned
                            // tree index forward to point at this new
                            // area of good beam

                            while (completeTime < macroTime-TIME_OFFSET)
                            {
                                do
                                {
                                    index++;
                                    tree->GetEntry(index);
                                    completeTime = (double)extTime*pow(2,32)+timetag;
                                }
                                while (chNo != detIndex);
                            }

                            cout << endl;

                            // now the cleaned tree index is pointing
                            // to the area of good beam, so continue
                            // filling the channel-specific tree
                        }

                        // lastly, because we shifted to a new
                        // macropulse at the start of this conditional,
                        // we should reset the event counter for this
                        // channel
                        evtNo = 0;
                    }

                    // now all DPP event variables have the correct values; fill the
                    // tree with an event
                    processDPPEvent(chTree);
                }

                else
                {
                    // we're in waveform mode - add an event to the waveform tree
                    processWaveformEvent(chTreeW);
                }

                // in preparation for looping to the next event, update counters
                evtNo++;
                prevCompleteTime = completeTime;
                prevEvtType = evtType;
            }
        }
    }
}

int main(int argc, char* argv[])
{

    // read in the raw tree name
    string runDir = argv[1];
    string runNo = argv[2];

    // Open the raw tree from the initial sort. If it doesn't exist, exit.
    TFile *file;

    stringstream treeName;
    stringstream fileName;

    treeName << runDir << "-" << runNo; 
    fileName << analysispath <<"analysis/" << runDir << "/" << treeName.str() << ".root";

    file = new TFile(fileName.str().c_str(),"UPDATE");

    if(file->Get("tree"))
    {
        // Found the raw tree; start sorting the raw trees into new trees for
        // each channel.
        cout << "Found previous raw tree; creating cleaned subtrees for " << treeName << endl;

        tree = (TTree*)file->Get("tree");

        // link the raw tree to our buffer variables
        tree->SetBranchAddress("chNo",&chNo);
        tree->SetBranchAddress("evtType",&evtType);
        tree->SetBranchAddress("timetag",&timetag);
        tree->SetBranchAddress("extTime",&extTime);
        tree->SetBranchAddress("sgQ",&sgQ);
        tree->SetBranchAddress("lgQ",&lgQ);
        tree->SetBranchAddress("fineTime",&fineTime);
        tree->SetBranchAddress("waveform",&dummyWaveform);

        ch0Tree = new TTree("ch0Tree","");
        branch(ch0Tree);

        ch2Tree = new TTree("ch2Tree","");
        branch(ch2Tree);

        ch4Tree = new TTree("ch4Tree","");
        branch(ch4Tree);

        ch6Tree = new TTree("ch6Tree","");
        branch(ch6Tree);

        ch0TreeW = new TTree("ch0TreeW","");
        branchW(ch0TreeW);

        ch2TreeW = new TTree("ch2TreeW","");
        branchW(ch2TreeW);

        ch4TreeW = new TTree("ch4TreeW","");
        branchW(ch4TreeW);
    }

    else
    {
        cout << "Failed to find raw tree - exiting " << fileName << endl;
        exit(1);
    }

    // We need to populate events into two trees for each channel, one for DPP
    // mode and one for waveform mode (except for ch6, which has no waveform
    // mode data).

    // populate target changer events
    processTargetChanger();

    // Now the macropulse structure is prepared in ch0Tree
    // link ch0Tree branches to macropulse variables in preparation for
    // filling ch2, ch4, ch6 trees
    ch0Tree->SetBranchAddress("macroNo",&macroNo);
    ch0Tree->SetBranchAddress("evtNo",&evtNo);
    ch0Tree->SetBranchAddress("macroTime",&macroTime);
    ch0Tree->SetBranchAddress("targetPos",&targetPos);

    // loop through raw tree multiple times, once per detector channel, and
    // populate detector events into channel-specific trees
    processDetEvents();

    file->Write();

    file->Close();
}
