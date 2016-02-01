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

// error log stream
ofstream error;


// EXPERIMENTAL PARAMETERS


// channel-specific time offset relative to the macropulse start time
double TIME_OFFSET;

// Target changer time delay after the macropulse start time as calculated using
// the summed detector signals' gamma peaks
const double MACROPULSE_OFFSET = 951; // in ns

// Time correction to MACROPULSE_OFFSET for the monitor; takes into account the
// additional ~10 meters between the summed detector and the extra cable length
const double MONITOR_OFFSET = 0; // in ns

// Time delay of scavenger afer summed det
const double SCAVENGER_OFFSET = 12.4; // in ns

// Period of macropulses
const double MACRO_PERIOD = 8330000; // in ns

// Target Changer lgQ gates for setting integral target positions 1-6
const int tarGate[12] = {5000,7500,11000,13500,17000,20000,23000,27000,29000,33000,35000,39000};



// TREE VARIABLES

// for holding event data from the raw tree
unsigned int chNo, evtType, extTime, fineTime, sgQ, lgQ;
double timetag;

// for holding event data for the new channel-specific trees
unsigned int evtNo, macroNo,/* microNo,*/ targetPos;
double completeTime, macroTime /*microTime*/;
vector<int> waveform;

// channel-specific trees' event structure (different from raw tree's)
struct event
{
  unsigned int macroNo; // label each event by macropulse
  //unsigned int microNo; // label each event by micropulse
  unsigned int evtNo; // uniquely label each event in a macropulse

  double macroTime; // the event's time-zero reference (the macropulse start)
  double completeTime; // the event's 48-bit timestamp
  //double microTime; // the event's timestamp relative to the micropulse it's in

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
    //tree->Branch("microNo",&ev.microNo,"microNo/i");
    tree->Branch("evtNo",&ev.evtNo,"evtNo/i");
    tree->Branch("macroTime",&ev.macroTime,"macroTime/d");
    tree->Branch("completeTime",&ev.completeTime,"completeTime/d");
    //tree->Branch("microTime",&ev.microTime,"microTime/d");
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
    //ev.microNo = microNo;
    ev.evtNo = evtNo;
    ev.macroTime = macroTime;
    ev.completeTime = completeTime;
    //ev.microTime = microTime;
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

    if (timeDiff > 0 && timeDiff < 650000)
    {
        getWaveform(dummyWaveform);
        fillTree(tree);
    }
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
    for(int k=0; k<6; k++)
    {
        if (lgQ>tarGate[2*k] && lgQ<tarGate[2*k+1])
        {
            // lgQ fits within one of the lgQ gates we defined, so
            // assign that target position
            return k+1;
        }
    }

    // lgQ doesn't fit within one of the lgQ gates, so set
    // the target position to 0. These will be thrown out
    // during plotting of cross-sections (no target was
    // in-beam during this event).
    return 0;
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

    // only use first 10% of entries (testing purposes)
    //totalEntries /= 10;

    cout << "Processing target changer events into ch0Tree." << endl;

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

            cout << "Target position = " << targetPos << ", macroNo = " << macroNo << ", macroTime = " << macroTime << "\r";
            fflush(stdout);

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

        // set the time offset to be applied to each channel based on which
        // channel we're sorting
        switch(detIndex)
        {
            case 2:
                chTree = ch2Tree;
                chTreeW = ch2TreeW;
                TIME_OFFSET = MACROPULSE_OFFSET+MONITOR_OFFSET;
                break;
            case 4:
                chTree = ch4Tree;
                chTreeW = ch4TreeW;
                TIME_OFFSET = MACROPULSE_OFFSET;
                break;
            case 6:
                chTree = ch6Tree;
                TIME_OFFSET = MACROPULSE_OFFSET+SCAVENGER_OFFSET;
                // no ch6 waveform tree
                break;
        }

        // all the preliminary steps ready:
        // start looping through detector events
        cout << "Processing channel " << detIndex << " events..." << endl;

        // use to loop through all entries and make sure we don't run off the
        // end of the raw tree
        int totalEntries = tree->GetEntries();

        // only take first 10% of events (testing mode)
        //totalEntries /= 10;

        // use to make sure we don't run off the end of the ch0Tree
        int ch0TreeEntries = ch0Tree->GetEntries();

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
                completeTime = (double)extTime*pow(2,32)+timetag+fineTime*(double)2/1024;
                //cout << completeTime << endl;

                // now check to see if the event is DPP or waveform mode
                if (evtType == 1)
                {
                    // DPP mode

                    // check to see if there was a timestamp-reset failure (the
                    // extended timestamp incremented, but the 32-bit timestamp
                    // didn't reset). If so, discard the event and move on.
                    if (completeTime > pow(2,33)-500 && completeTime < pow(2,33)+500)
                    {
                        error << "Found an event with a timestamp-reset failure (i.e., extTime incremented before timetag was reset to 0). Macrotime = " << macroTime << ", macroNo = " << macroNo << "completeTime = " << completeTime << endl;
                        error << "Skipping to next event..." << endl;
                        continue;
                    }

                    // Now check to see if we've wrapped around to a new DPP 
                    // mode period.
                    // When a DPP/waveform mode switch occurs, the timestamps are reset to 0 at
                    // the start of the new mode. To keep track of when a channel switches
                    // modes, we need to compare consecutive events' timestamps.
                    if (completeTime < prevCompleteTime)
                    {
                        // new DPP mode period, triggered by a
                        // channel-specific event wrapping back around
                        // to time zero

                        // update the ch0Tree to point at the new
                        // macropulse
                        do
                        {
                            if (macroNo+1>=ch0TreeEntries)
                            {
                                cout << "Exceeded the last macropulse on ch0Tree. Exiting loop." << endl;
                                break;
                            }

                            prevMacroTime = macroTime;
                            ch0Tree->GetEntry(macroNo+1);
                        }
                        while (prevMacroTime < macroTime);

                        if (macroNo+1>=ch0TreeEntries)
                        {
                            cout << "Exceeded the last macropulse on ch0Tree. Exiting loop." << endl;
                            break;
                        }

                        cout << "Start new DPP period @ macroNo " << macroNo << "\r";
                        fflush(stdout);
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

                        if (macroNo+1>=ch0TreeEntries)
                        {
                            cout << "Exceeded the last macropulse on ch0Tree. Exiting loop." << endl;
                            break;
                        }
                        // we're not at the end of the run, so we
                        // should pull the next macropulse to see
                        // what to do with it

                        ch0Tree->GetEntry(macroNo+1);

                        //cout << "macropulse expired; new macroTime = " << macroTime << ", macroNo = " << macroNo << endl;
                        //fflush(stdout);

                        // evtNo was updated when we pulled the next
                        // macropulse - let's check to see if we've
                        // reached the end of a DPP period
                        if (evtNo==1)
                        {
                            // new DPP period found because the
                            // macropulse has evtNo==1, so move the
                            // cleaned tree indices forward to the new
                            // DPP period

                            do
                            {
                                do
                                {
                                    if (index==totalEntries)
                                    {
                                        cout << "Reached the end of the raw tree; end looping." << endl;
                                        break;
                                    }
                                    index++;
                                    tree->GetEntry(index);
                                }
                                while (chNo != detIndex || evtType != 1);

                                if (index==totalEntries)
                                    {
                                        break;
                                    }

                                prevCompleteTime = completeTime;
                                completeTime = (double)extTime*pow(2,32)+timetag;
                            }
                            while (prevCompleteTime < completeTime);

                            if (index==totalEntries)
                            {
                                break;
                            }

                            cout << "Start new DPP period @ macroNo " << macroNo << "\r";
                            fflush(stdout);
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

                            // keep track of the first macroTime with a beam
                            // anomaly so we can see if the anomaly spans a mode
                            // change
                            double beamAnomalyStart = macroTime;

                            error << "Macropulse timing anomaly found at macroTime = " << beamAnomalyStart << endl;
                            error << "prevMacroTime = " << ", macroNo = " << macroNo << ", completeTime = " << completeTime << endl;

                            while(!(macroTime-prevMacroTime > 8300000 && macroTime-prevMacroTime < 8360000) && !(macroTime-prevMacroTime > 16600000 && macroTime-prevMacroTime < 17200000) && !(macroTime-prevMacroTime > 24900000 && macroTime-prevMacroTime < 25080000))
                            {
                                // update the ch0Tree to point at the
                                // next macropulse so we can check to
                                // see if it's out-of-sync with the
                                // previous macropulse 

                                //cout << "macroNo = " << macroNo << ", macroTime = " << macroTime << ", prevMacroTime = " << prevMacroTime << endl;

                                if (macroNo+1>=ch0TreeEntries)
                                {
                                    cout << "Exceeded the last macropulse on ch0Tree. Exiting loop." << endl;
                                    break;
                                }

                                prevMacroTime = macroTime;
                                ch0Tree->GetEntry(macroNo+1);

                            };

                            error << "Beam anomaly finished; new macroTime is " << macroTime << ", macroNo " << macroNo << ", completeTime = " << completeTime << endl;

                            // OK - ch0Tree is now in a period of
                            // good beam. We need to move the cleaned
                            // tree index forward to point at this new
                            // area of good beam

                            // Two steps must be taken to make sure
                            // we've moved the channel index far enough
                            // forward:
                            //
                            // 1) if the beam anomaly spanned a
                            // switch from DPP to waveform mode, we
                            // must first move the cleaned tree index forward to
                            // the new DPP period
                            //
                            // 2) move the cleaned tree index forward until the
                            // detector events are after the time-zero of the
                            // current macropulse

                            // test for step 1
                            if (beamAnomalyStart > macroTime)
                            {
                                // beam anomaly DID span a DPP/waveform mode
                                // change
                                error << "Beam anomaly spanned DPP/waveform mode change; moving the channel index forward to the next DPP mode period." << endl;
                                do
                                {
                                    do
                                    {
                                        if (index==totalEntries)
                                        {
                                            error << "Reached the end of the raw tree; end looping." << endl;
                                            break;
                                        }

                                        index++;
                                        tree->GetEntry(index);
                                    }
                                    while (chNo != detIndex || evtType != 1);

                                    if (index==totalEntries)
                                    {
                                        break;
                                    }

                                    prevCompleteTime = completeTime;
                                    completeTime = (double)extTime*pow(2,32)+timetag;
                                }
                                // keep looping until the timestamps reset to 0 for
                                // this channel (i.e., new DPP period)
                                while (prevCompleteTime < completeTime);
                            }

                            error << "Finished moving channel index forward to the next DPP mode period." << endl;

                            // if we've reached the end of the raw tree events,
                            // exit the loop
                            if (index==totalEntries)
                            {
                                break;
                            }

                            // continue moving index forward until we reach the
                            // current macropulse time
                            do
                            {
                                do
                                {
                                    index++;
                                    tree->GetEntry(index);
                                }
                                while (chNo != detIndex || evtType != 1);

                                completeTime = (double)extTime*pow(2,32)+timetag;
                            }
                            while (completeTime < macroTime-TIME_OFFSET);

                            // now the cleaned tree index is pointing
                            // to the area of good beam, so continue
                            // filling the channel-specific tree
                        }

                        // because we shifted to a new
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

    cout.precision(13);

    // read in the raw tree name
    string runDir = argv[1];
    string runNo = argv[2];

    // Open the raw tree from the initial sort. If it doesn't exist, exit.
    TFile *file, *fileOut;

    stringstream treeName;
    stringstream fileInName;
    stringstream fileOutName;
    stringstream errorName;

    treeName << runDir << "-" << runNo; 

    // write out sorting errors to errors.log
    errorName << analysispath <<  "analysis/" << runDir << "/" << treeName.str() << "_error.log";
    error.open(errorName.str());

    fileInName << analysispath <<"analysis/" << runDir << "/" << treeName.str() << "_raw.root";
    fileOutName << analysispath <<"analysis/" << runDir << "/" << treeName.str() << "_sorted.root";

    fileOut = new TFile(fileOutName.str().c_str(),"UPDATE");

    if(fileOut->Get("ch0Tree"))
    {
        cout << "Found previously existing channel-specific sorted trees - skipping channel-specific sorting." << endl;
        exit(1);
    }

    else
    {
        // no pre-existing output file: we'll need to do a new sort on the raw
        // tree data. So, first create the channel-specific trees:
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

        file = new TFile(fileInName.str().c_str(),"READ");

        if(file->Get("tree"))
        {
            // Found the raw tree; start sorting the raw trees into new trees for
            // each channel.
            cout << "Found a raw event tree; creating cleaned subtrees." << endl;

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
        }

        else
        {
            cout << "Failed to find raw tree - exiting " << fileInName << endl;
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

        fileOut->Write();

        fileOut->Close();
    }
}
