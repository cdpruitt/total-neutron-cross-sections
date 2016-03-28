/******************************************************************************
                                    RESORT 
******************************************************************************/
// This sub-code takes a raw ROOT tree (from ./raw) event data and processes it
// into channel-specific trees useful for analysis.

// Once the raw tree is opened, target changer events are passed into a
// target-changer-specific sub-tree and are labeled by macropulse number.
// Then, we loop through the raw tree again, pulling out events in channels 2
// (the monitor) and channel 4 (the detector). These events' timestamps are
// compared to the timestamps in the target changer sub-tree, allowing us to
// assign events in channels 2 and 4 to the appropriate macropulse.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TEntryList.h"

using namespace std;

/******************************************************************************/
// DEFINE EXPERIMENTAL CONSTANTS

// channel-specific time offset relative to the macropulse start time
double TIME_OFFSET;

// Target changer time delay after the macropulse start time; determined by
// adjusting until the gamma peaks appear in the right place for the given
// target distance from the neutron source
const double MACROPULSE_OFFSET = 842; // in ns
// 814 for Ne 

// Time correction to MACROPULSE_OFFSET for the monitor; takes into account the
// additional ~10 meters between the summed detector and the extra cable length
// NOTE: currently unused - uncomment to use
const double MONITOR_OFFSET = 0; // in ns

// Time delay of scavenger afer summed det
const double SCAVENGER_OFFSET = -12.51; // in ns

// Beam macropulse period
const double MACRO_PERIOD = 8333000; // in ns

// Target Changer lgQ gates for setting integral target positions 1-6
const int tarGate[12] = {5000,7500,11000,13500,17000,20000,23000,27000,29000,33000,35000,39000};


/******************************************************************************/
// EVENT VARIABLES

// for holding event data from the raw tree
unsigned int chNo, evtType, extTime, fineTime, sgQ, lgQ;
double timetag;

// for holding event data for the new channel-specific trees
unsigned int evtNo, macroNo, targetPos;
double completeTime, macroTime;
vector<int> waveform;

// channel-specific trees' event structure (different from raw tree's)
struct event
{
  unsigned int macroNo; // label each event by macropulse
  unsigned int evtNo; // uniquely label each event in a macropulse

  double macroTime; // the event's time-zero reference (the macropulse start)
  double completeTime; // the event's 48-bit timestamp

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
double extTimePrev = 0;
double timetagPrev = 0;

int prevEvtType;

// detector index (2 is monitor, 4 is summed detector, 6 is scavenger)
int detIndex;

/******************************************************************************/
// TREES

// Create pointer to the already-populated TTree 'tree' from ./raw
// We will loop through this tree and divide up its contents by channel to
// perform later analysis (like assigning TOFs, calculating cross-sections, etc)
TTree *tree;

// Create pointers to the to-be-created channel-specific sub-trees.
// Each channel has a tree for DPP and waveform mode
TTree* ch0Tree;
TTree* ch2Tree;
TTree* ch4Tree;
TTree* ch6Tree;
TTree* ch0TreeW;
TTree* ch2TreeW;
TTree* ch4TreeW;


/******************************************************************************/
// DEBUGGING

// If errors are encountered during the resort, print them out to an error file
ofstream error;


/******************************************************************************/
// HELPER METHODS

// Used to connect a channel-specific tree to DPP event variables so we can
// start populating it with DPP events
void branch(TTree* tree)
{
    tree->Branch("macroNo",&ev.macroNo,"macroNo/i");
    tree->Branch("evtNo",&ev.evtNo,"evtNo/i");
    tree->Branch("macroTime",&ev.macroTime,"macroTime/d");
    tree->Branch("completeTime",&ev.completeTime,"completeTime/d");
    tree->Branch("targetPos",&ev.targetPos,"targetPos/i");
    tree->Branch("sgQ",&ev.sgQ,"sgQ/i");
    tree->Branch("lgQ",&ev.lgQ,"lgQ/i");
    tree->Branch("waveform",&ev.waveform);
}

// Used to connect a channel-specific tree to waveform event variables so we can
// start populating it with waveform events
void branchW(TTree* tree)
{
    tree->Branch("macroNo",&ev.macroNo,"macroNo/i");
    tree->Branch("evtNo",&ev.evtNo,"evtNo/i");
    tree->Branch("completeTime",&ev.completeTime,"completeTime/d");
    tree->Branch("targetPos",&ev.targetPos,"targetPos/i");
    tree->Branch("waveform",&ev.waveform);
}

// add an event to a channel-specific tree. This is called in the main event
// loop once all the event variables have been filled.
void fillTree(TTree* tree)
{
    ev.macroNo = macroNo;
    ev.evtNo = evtNo;
    ev.macroTime = macroTime;
    ev.completeTime = completeTime;
    ev.targetPos = targetPos;
    ev.sgQ = sgQ;
    ev.lgQ = lgQ;

    // To save only 1 DPP event's waveform in every 10000 events, uncomment this
    // (Note - will reduce output tree size by ~half)
    /*if ((int)completeTime%10000 != 0 && evtType==1)
    {
        waveform.clear(); 
    }*/

    ev.waveform = waveform; 

    tree->Fill();
}

// access an event's waveform data in preparation for storing it in a
// channel-specific tree
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
    double timeDiff = completeTime-macroTime;

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

    // DPP mode is event type 1; waveform mode is event type 2
    // We need to check and see if this waveform event is the first after a
    // mode switch
    if (prevEvtType==1)
    {
        // A mode switch has happened: reset the evtNo so that waveforms in
        // this macropulse will be labeled 0, 1, 2, etc.
        evtNo=0;
    }

    fillTree(tree);
}


// Use the lgQ from the target changer to assign each macropulse a target
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

/******************************************************************************/
// MAIN METHODS (i.e., used to directly process events)


// fill the ch0Tree with target changer events so that we have macropulse
// start times and target positions available to assign to the detector
// channels
// ch0Tree is the first channel-specific tree we fill
void processTargetChanger()
{
    // Populate a tree of only target changer events in preparation
    // for assigning times and target changer positions to the
    // detector channels (ch2, ch4, ch6)

    // macroNo numbering starts at 0 (start of the sub-run)
    macroNo = 0;

    // loop through raw tree to find target changer events
    int totalEntries = tree->GetEntries();

    // only use first 10% of entries (for debugging)
    // uncomment to use
    //totalEntries /= 10;

    cout << "Processing target changer events into ch0Tree." << endl;

    // Main loop through the events of the raw tree
    for (int i=0; i<totalEntries; i++)
    {
        tree->GetEntry(i);

        if (chNo == 0)
        {
            if (lgQ == 65535)
            {
                // ignore target changer events that with off-scale integrated
                // charge - these are suspected to be retriggers and are NOT
                // time-correlated with the facility RF time reference.
                continue;
            }

            // assign the macropulse time reference (start time) using the
            // target changer event timestamps
            macroTime = (double)extTime*pow(2,32)+timetag;

            // assign this event's timestamp (because this is a target changer
            // event, this is the same as the macroTime)
            completeTime = macroTime;

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
                    // NOT the start of a new DPP period (99.9% of the time)
                    evtNo = 0;
                }

                if (extTime > extTimePrev && timetag > pow(2,32)-10000)
                {
                    error << "Found a target changer event with a timestamp-reset failure (i.e., extTime incremented before timetag was reset to 0). MacroNo = " << macroNo << ", completeTime = " << completeTime << ", extTime = " << extTime << ", extTimePrev = " << extTimePrev << ", timetag = " << timetag << ", timetagPrev = " << timetagPrev << endl;
                    error << "Skipping to next target changer event..." << endl;
                    continue;
                }

                // figure out which target is in the beam during this event
                targetPos = assignTargetPos(lgQ);

                // all the variables are updated, so fill the target changer
                // tree with the event
                fillTree(ch0Tree);

                // increment macropulse counter to prep for next macropulse
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

            // before we pull the next event, save the previous event's data so we can
            // compare them with the new event's data
            prevEvtType = evtType;
            extTimePrev = extTime;
            timetagPrev = timetag;

            cout << "Target position = " << targetPos << ", macroNo = " << macroNo << ", macroTime = " << macroTime << "\r";
            fflush(stdout);
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

        // each event should be uniquely labeled by its channel-specific order
        // in the macropulse
        // reset event counter in preparation for looping through events
        evtNo = 0; 

        // point to the first macropulse in preparation for assigning events to
        // macropulses
        ch0Tree->GetEntry(0);

        // in preparation for comparing reset the timestamp holders
        prevMacroTime = 0;
        prevCompleteTime = 0;
        extTimePrev = 0;
        timetagPrev = 0;

        // Define 'holders' for the channel-specific trees so that we can
        // abstract tree processing onto the holders rather then directly on
        // the trees themselves
        TTree *chTree, *chTreeW;

        // Each channel has a different time offset; set the time offset we
        // should apply, depending on channel
        // Also, assign the correct tree to the 'holders' just defined
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
            default:
                cout << "Error: non-existent channel index given by detIndex" << endl;
                exit(1);
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

        // Main loop through monitor/detector events
        for (int index = 0; index<totalEntries; index++)
        {
            // pull next event's data from the raw tree
            tree->GetEntry(index);

            // Only process events from the correct channel:
            if (chNo == detIndex)
            {
                // we've located a detector event from the correct (detIndex) channel.

                // use its data along with the target changer data to assign this
                // event the correct time, macropulse number, etc.

                // first calculate this event's time
                completeTime = (double)extTime*pow(2,32)+timetag+fineTime*(double)2/1024+TIME_OFFSET;

                // now check to see if the event is DPP or waveform mode
                if (evtType == 1)
                {
                    // DPP mode

                    // Now, it's time to apply some filters to ignore events
                    // that have problems with their time structure

                    // Check to see if there was a timestamp-reset failure (the
                    // extended timestamp incremented, but the 32-bit timestamp
                    // didn't reset). If so, discard the event and move on.

                    if (extTime > extTimePrev && timetag > pow(2,32)-10000)
                    {
                        error << endl << "MACRO NO " << macroNo << endl;
                        error << "Found a detector event with a timestamp-reset failure (i.e., extTime incremented before timetag was reset to 0). completeTime = " << completeTime << endl;
                        error << "Skipping to next event..." << endl;
                        continue;
                    }

                    // Next, check to see if this is event is the first one
                    // of a new DPP mode period.

                    // When a DPP/waveform mode switch occurs, the timestamps are reset to 0 at
                    // the start of the new mode. To keep track of when a channel switches
                    // modes, we should compare consecutive events' timestamps.
                    if (completeTime < prevCompleteTime)
                    {
                        // new DPP mode period, triggered by this event's
                        // timestamp wrapping back around to time zero

                        // move to the next macropulse of the ch0Tree
                        do
                        {
                            if (macroNo+1>=ch0TreeEntries)
                            {
                                // ran out of new macropulses - start a series
                                // of breaks to exit the main loop
                                cout << "Exceeded the last macropulse on ch0Tree. Exiting loop." << endl;
                                break;
                            }

                            prevMacroTime = macroTime;
                            ch0Tree->GetEntry(macroNo+1);
                        }
                        while (prevMacroTime < macroTime);

                        if (macroNo+1>=ch0TreeEntries)
                        {
                            // next break statement, needed to exit the main loop
                            break;
                        }

                        // ch0Tree is now pointing at the first
                        // macropulse of the next DPP mode period
                        cout << "Start new DPP period @ macroNo " << macroNo << "\r";
                        fflush(stdout);
                        
                        // reset the event counter for this channel
                        evtNo = 0;
                    }

                    // check to see if enough time has elapsed that the current
                    // macropulse has expired (i.e., beam has gone off again)
                    else if (completeTime-macroTime-TIME_OFFSET > 8000000)
                    {
                        // macropulse expired
                        // time to switch to the next macropulse

                        // hold the old macro time so we'll be able to
                        // compare it with the new macro time and make
                        // sure there's still normal beam structure
                        // (i.e., no beam-off periods or accidental
                        // re-triggerings)
                        prevMacroTime = macroTime;

                        if (macroNo+1>=ch0TreeEntries)
                        {
                            cout << "Exceeded the last macropulse on ch0Tree. Exiting loop." << endl;
                            break;
                        }

                        // we're not at the end of the run, so proceed to the
                        // next macropulse on ch0Tree
                        ch0Tree->GetEntry(macroNo+1);

                        // If the macropulse was the first macropulse of a new
                        // DPP mode, then evtNo will have been updated
                        // to be 1 (otherwise, it's 0).
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
                                completeTime = (double)extTime*pow(2,32)+timetag+TIME_OFFSET;
                            }
                            while (prevCompleteTime < completeTime || extTime > 0);

                            if (index==totalEntries)
                            {
                                break;
                            }

                            cout << "Start new DPP period @ macroNo " << macroNo << "\r";
                            fflush(stdout);
                        }

                        // not a new DPP period
                        // Lastly, let's check to see if there are any beam
                        // anomalies (that is, the new macropulse is
                        // out-of-sync with the RF frequency )
                        else if (!(macroTime-prevMacroTime > MACRO_PERIOD-30000 && macroTime-prevMacroTime < MACRO_PERIOD+30000) && !(macroTime-prevMacroTime > 2*(MACRO_PERIOD-30000) && macroTime-prevMacroTime < 2*(MACRO_PERIOD+30000)) && !(macroTime-prevMacroTime > 3*(MACRO_PERIOD-30000) && macroTime-prevMacroTime < 3*(MACRO_PERIOD+30000)) && !(macroTime-prevMacroTime > 4*(MACRO_PERIOD-30000) && macroTime-prevMacroTime < 4*(MACRO_PERIOD+30000)))
                        {
                            // there is a beam anomaly - cycle through
                            // the ch0Tree until we find
                            // another period of clean beam

                            // keep track of the first macroTime with a beam
                            // anomaly so we can see if the anomaly spans a mode
                            // change
                            double beamAnomalyStart = macroTime;

                            error << endl << "MACRO NO " << macroNo << endl;
                            error << "Macropulse timing anomaly found at macroTime = " << beamAnomalyStart << endl;
                            error << "macroTime - prevMacroTime = " << macroTime-prevMacroTime << endl;
                            error << "Looping through macropulse list..." << endl;

                            do
                            {
                                // update the ch0Tree to point at the
                                // next macropulse so we can check to
                                // see if it's out-of-sync with the
                                // previous macropulse 

                                if (macroNo+1>=ch0TreeEntries)
                                {
                                    cout << "Exceeded the last macropulse on ch0Tree. Exiting loop." << endl;
                                    break;
                                }

                                prevMacroTime = macroTime;
                                ch0Tree->GetEntry(macroNo+1);
                                error << "macroTime = " << macroTime << "; macroTime-prevMacroTime = " << macroTime-prevMacroTime << endl;
                            }
                            while (!(macroTime-prevMacroTime > MACRO_PERIOD-30000 && macroTime-prevMacroTime < MACRO_PERIOD+30000) && !(macroTime-prevMacroTime > 2*(MACRO_PERIOD-30000) && macroTime-prevMacroTime < 2*(MACRO_PERIOD+30000)) && !(macroTime-prevMacroTime > 3*(MACRO_PERIOD-30000) && macroTime-prevMacroTime < 3*(MACRO_PERIOD+30000)) && !(macroTime-prevMacroTime > 4*(MACRO_PERIOD-30000) && macroTime-prevMacroTime < 4*(MACRO_PERIOD+30000)));

                            error << "Beam anomaly finished; new macroTime is " << macroTime << ", macroNo " << macroNo << ", completeTime = " << completeTime << endl;

                            // OK - ch0Tree is now pointing at a macropulse in a
                            // period of good beam. We need to move forward in
                            // the raw tree to point at the same area of good
                            // beam

                            // Two steps must be taken to make sure
                            // we've moved the channel index far enough
                            // forward:
                            //
                            // 1) if the beam anomaly spanned a
                            // switch from DPP to waveform mode, we
                            // must first move the raw tree index forward to
                            // the new DPP mode's events
                            //
                            // 2) move the raw tree index forward until the
                            // detector events are after the time-zero of the
                            // current macropulse

                            // test for step 1
                            if (beamAnomalyStart > macroTime)
                            {
                                // beam anomaly DID span a DPP/waveform mode
                                // move the raw tree forward past mode change
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
                                        // Reached the end of raw tree;
                                        // end sorting.
                                        break;
                                    }

                                    prevCompleteTime = completeTime;
                                    completeTime = (double)extTime*pow(2,32)+timetag+TIME_OFFSET;
                                }

                                // keep looping until the timestamps reset to 0 for
                                // this channel (i.e., new DPP period)
                                while (prevCompleteTime < completeTime);

                                error << "Finished moving channel index forward to the next DPP mode period." << endl;
                            }

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

                                completeTime = (double)extTime*pow(2,32)+timetag+TIME_OFFSET;
                            }
                            while (completeTime < macroTime);

                            // now the raw tree index is pointing
                            // to the area of good beam, so continue
                            // filling the channel-specific tree
                        }

                        // because we shifted to a new
                        // macropulse at the start of this conditional,
                        // we should reset the event counter for this
                        // channel
                        evtNo = 0;
                    }

                    // now all DPP event variables have the correct values, and
                    // we're pointing at the correct macropulse in ch0Tree
                    // --> fill the channel-specific tree with the event
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
                extTimePrev = extTime;
                timetagPrev = timetag;
            }
        }
    }
}

/******************************************************************************/
// MAIN

// Set cout precision to display the entire value of event timestamps
int main(int argc, char* argv[])
{
    cout << endl << "Entering re-sort..." << endl;

    cout.precision(13);
    error.precision(13);

    // read in the raw tree name
    string runDir = argv[1];
    string runNo = argv[2];
    string outpath = argv[3];

    // Open the raw tree from the initial sort. If it doesn't exist, exit.
    TFile *file, *fileOut;

    stringstream treeName;
    stringstream fileInName;
    stringstream fileOutName;
    stringstream errorName;

    treeName << "run" << runDir << "-" << runNo; 

    // Create an error log where sorting errors are recorded
    errorName << outpath <<  "/analysis/run" << runDir << "/" << treeName.str() << "_error.log";
    error.open(errorName.str());

    // Prepare filenames for input file (raw tree) and output file (channel-
    // specific trees)
    fileInName << outpath <<"/analysis/run" << runDir << "/" << treeName.str() << "_raw.root";
    fileOutName << outpath <<"/analysis/run" << runDir << "/" << treeName.str() << "_sorted.root";

    // Check to see if an output file already exists. If it does, abort sorting
    // and pass control back up to ./sort.sh
    fileOut = new TFile(fileOutName.str().c_str(),"UPDATE");
    if(fileOut->Get("ch0Tree"))
    {
        cout << fileOutName.str() << " already exists. Skipping ./resort..." << endl;
        exit(0);
    }

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

    if(!file->Get("tree"))
    {
        cout << "Error: failed to find raw tree in " << fileInName.str() << endl;
        exit(1);
    }

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
