/******************************************************************************
                                    RESORT 
******************************************************************************/
// This sub-code takes a ROOT tree of raw events (produced by running ./raw),
// splits it into channel-specific ROOT trees, and assigns each event to a
// macropulse in preparation for producing cross section plots.

// This process involves two steps:
// 1. Split target changer events into a separate target-changer-specific tree
// where they are labeled by macropulse number.
// 2. Using the target-changer tree to keep track of which macropulse we're in,
// assign macropulses to channel 2 (the flux monitor) and channel 4 (the
// detector).

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

// Time offsets:

// "Macropulse offset" is the time difference between the target changer time
// and the true the macropulse start time. It is determined by trial-and-error
// until it moves the gamma peaks to the right place for the given
// target distance from the neutron source
const double MACROPULSE_OFFSET = 842; // in ns
// 814 for Ne 

// "Scavenger offset" is the time delay of scavenger relative to the main
// detector channel
const double SCAVENGER_OFFSET = -12.51; // in ns

// "Macropulse period" is the time between macropulses (currently 120 Hz at WNR)
const double MACRO_PERIOD = 8333300; // in ns

// "Macropulse length" is (micropulse period)*(# micropulses per macropulse)
const double MACROPULSE_LENGTH = 650000; // in ns

// "Sync window" defines the acceptable time variation between consecutive macropulses
// (with respect to the macropulse period) in order for those macropulses to be
// considered still "in sync" with the macropulse period
const double SYNC_WINDOW = 0.002; // as fraction of MACRO_PERIOD 

// "Target changer charge gates" are used to assign the target changer position
// based on the target changer signal's integrated charge
const int tarGate[12] = {5000,7500,11000,13500,17000,20000,23000,27000,29000,33000,35000,39000};
            // Position: 1low  1hi  2low   2hi  3low   3hi  4low   4hi  5low   5hi  6low   6hi

// Establish which channels are active in each mode
const vector<string> activeDPPChannels = {"ch0","ch2","ch4","ch6"};
const vector<string> activeWaveformChannels = {"ch0","ch2","ch4"};

const double SCALEDOWN = 1; // (for debugging) only sort (total/SCALEDOWN) events

/******************************************************************************/
// track number of processed events of each type
long numberOfDPPs = 0;
long numberOfWaveforms = 0;
long numberOfCh0Waveforms = 0;
long numberOfCh2Waveforms = 0;
long numberOfCh4Waveforms = 0;

/******************************************************************************/
// Define tree structures and functions to fill/read from trees
/******************************************************************************/

// used to sort events into channel-specific trees, but do no processing
struct RawEvent
{
    double timetag; // 1 sample granularity, 32 bits
    unsigned int extTime;
    unsigned int fineTime; // provide additional bits of granularity 
    unsigned int evtNo;
    unsigned int evtType;
    unsigned int chNo;

    unsigned int sgQ; // integrated charge of event from short gate
    unsigned int lgQ; // integrated charge of event from long gate
    vector<int> *waveform; // contains all waveform samples for each event to allow for corrections in analysis
} rawEvent;

// used to store processed events after they have been mated with a macropulse
struct ProcessedEvent
{
    unsigned int macroNo; // label each event by macropulse
    double macroTime; // provide the macropulse "zero" time
    unsigned int evtNo; // uniquely label each event in a macropulse
    double completeTime; // the event's 48-bit timestamp
    unsigned int targetPos; // target position
    unsigned int sgQ, lgQ; // the event's short and long integrated charge gates
    vector<int> *waveform; // waveform data for this event
} procEvent;

// used to keep track of the macropulse structure of the sub-run.
struct TargetChangerEvent
{
    unsigned int macroNo; // label each event by macropulse
    unsigned int modeChange; // indicate the first event after a mode change
    double macroTime; // the event's time-zero reference (the macropulse start)
    unsigned int targetPos; // target position
} tcEvent;

// Used to connect a channel-specific tree to DPP event variables so we can
// start populating it with DPP events
void branchRaw(TTree* tree)
{
    tree->Branch("timetag",&rawEvent.timetag,"timetag/d");
    tree->Branch("extTime",&rawEvent.extTime,"extTime/i");
    tree->Branch("fineTime",&rawEvent.fineTime,"fineTime/i");
    tree->Branch("sgQ",&rawEvent.sgQ,"sgQ/i");
    tree->Branch("lgQ",&rawEvent.lgQ,"lgQ/i");
    tree->Branch("waveform",&rawEvent.waveform);
}

// Used to connect a channel-specific tree to waveform event variables so we can
// start populating it with waveform events
void branchRawW(TTree* tree)
{
    tree->Branch("timetag",&rawEvent.timetag,"timetag/d");
    tree->Branch("extTime",&rawEvent.extTime,"extTime/i");
    tree->Branch("evtNo",&rawEvent.evtNo,"evtNo/i");
    tree->Branch("waveform",&rawEvent.waveform);
}

// Used to connect a channel-specific tree to DPP event variables so we can
// start populating it with DPP events
void branchProc(TTree* tree)
{
    tree->Branch("macroNo",&procEvent.macroNo,"macroNo/i");
    tree->Branch("macroTime",&procEvent.macroTime,"macroTime/d");
    tree->Branch("evtNo",&procEvent.evtNo,"evtNo/i");
    tree->Branch("completeTime",&procEvent.completeTime,"completeTime/d");
    tree->Branch("targetPos",&procEvent.targetPos,"targetPos/i");
    tree->Branch("sgQ",&procEvent.sgQ,"sgQ/i");
    tree->Branch("lgQ",&procEvent.lgQ,"lgQ/i");
    tree->Branch("waveform",&procEvent.waveform);
}

// Used to connect a channel-specific tree to waveform procEvent variables so we can
// start populating it with waveform procEvents
void branchProcW(TTree* tree)
{
    tree->Branch("macroNo",&procEvent.macroNo,"macroNo/i");
    tree->Branch("evtNo",&procEvent.evtNo,"evtNo/i");
    tree->Branch("completeTime",&procEvent.completeTime,"completeTime/d");
    tree->Branch("targetPos",&procEvent.targetPos,"targetPos/i");
    tree->Branch("waveform",&procEvent.waveform);
}

void branchTargetChanger(TTree* tree)
{
    tree->Branch("macroNo",&tcEvent.macroNo,"macroNo/i");
    tree->Branch("macroTime",&tcEvent.macroTime,"macroTime/d");
    tree->Branch("modeChange",&tcEvent.modeChange,"modeChange/i");
    tree->Branch("targetPos",&tcEvent.targetPos,"targetPos/i");
}

void fillRawTree(TTree* tree)
{
    // To save only 1 DPP event's waveform in every 10000 events, uncomment this
    // (Note - will reduce output tree size by ~half)
    /*if ((int)completeTime%10000 != 0 && evtType==1)
      {
      waveform->clear(); 
      }*/
    tree->Fill();

    if(rawEvent.evtType==2)
    {
        if(rawEvent.chNo==0)
        {
            numberOfCh0Waveforms++;
        }

        if(rawEvent.chNo==2)
        {
            numberOfCh2Waveforms++;
        }

        if(rawEvent.chNo==4)
        {
            numberOfCh4Waveforms++;
        }
    }
}

// add an event to a channel-specific tree. This is called in the main event
// loop once all the event variables have been filled.
void fillProcessedTree(TTree* tree)
{
    // To save only 1 DPP event's waveform in every 10000 events, uncomment this
    // (Note - will reduce output tree size by ~half)
    /*if ((int)completeTime%10000 != 0 && evtType==1)
    {
        waveform->clear(); 
    }*/
    tree->Fill();
}
/******************************************************************************/



/******************************************************************************/
// Event processing functions
/******************************************************************************/

// access an event's waveform data in preparation for storing it in a
// channel-specific tree
void getWaveform(vector<int>* dummyWaveform)
{
    // empty the stale waveform data
    rawEvent.waveform->clear();

    // fill waveform with wavelet data from an event
    // (the raw tree's branch points to dummyWaveform)
    for(int k=0; (size_t)k<dummyWaveform->size(); k++)
    {
        rawEvent.waveform->push_back(dummyWaveform->at(k));
    }
}

// Use the lgQ from the target changer to determine the target position
int assignTargetPos(int lgQ)
{
    for(int k=0; k<6; k++)
    {
        if (lgQ>tarGate[2*k] && lgQ<tarGate[2*k+1])
        {
            // lgQ fits within one of the lgQ gates we defined
            return k+1;
        }
    }

    // lgQ doesn't fit within one of the lgQ gates, so set
    // the target position to 0 (discarded in later analysis).
    return 0;
}

// Populate events from the input tree into channel-specific trees.
void separateByChannel(TTree* inputTree, vector<TTree*> orchardRaw, vector<TTree*> orchardRawW)
{
    cout << "Separating events by channel and event type..." << endl;

    int prevEvtType = 0;
    
    int totalEntries = inputTree->GetEntries();

    // Uncomment to sort only first 10% of events (debugging mode)
    totalEntries /= SCALEDOWN;

    // To uniquely identify each event, we assign each channel's events an event
    // number (evtNo), which is the event's order in its macropulse.
    int evtNo[activeDPPChannels.size()];
    memset(evtNo,0,sizeof(evtNo));
    
    // Loop through all events in input tree and separate them into channel-
    // specific trees
    for (int index=0; index<totalEntries; index++)
    {
        inputTree->GetEntry(index);

        // if the event is a target changer event (chNo = 0),
        // processTargetChanger will handle it separately because target
        // changer events are special, non-detector events
        if(rawEvent.chNo==0)
        {
            continue;
        }

        if (rawEvent.evtType==1)
        {
            // DPP mode
            fillRawTree(orchardRaw[rawEvent.chNo/2]);
        }

        if (rawEvent.evtType==2)
        {
            // Waveform mode
            // Because waveform mode events are not time referenced to target
            // changer events, we need to manually increment the evtNo here 
            if (prevEvtType==1)
            {
                // New waveform mode period
                memset(evtNo,0,sizeof(evtNo));
            }

            rawEvent.evtNo = evtNo[rawEvent.chNo/2];
            evtNo[rawEvent.chNo/2]++;

            fillRawTree(orchardRawW[rawEvent.chNo/2]);
        }

        if(index%10000==0)
        {
            cout << "processed " << index << " events.\r";
            fflush(stdout);
        }

        prevEvtType = rawEvent.evtType;
    }

    cout << "Separated " << totalEntries << " events into channels 2, 4 and 6." << endl;
}

// assign a macropulse, target position, etc to every target changer event
void processTargetChanger(TTree* inputTree, TTree* targetChangerTree, ofstream& error)
{
    // Define dummy variables for comparing the timestamps of consecutive events.
    // This will allow us to track time resets and DPP/waveform mode changes.
    double extTimePrev = 0;
    double timetagPrev = 0;
    int prevEvtType = 0;

    cout << "Processing target changer events into targetChangerTree." << endl;

    // The first macropulse is defined as 0 (start of the sub-run)
    tcEvent.macroNo = 0;

    // loop through raw tree to find target changer events
    int totalEntries = inputTree->GetEntries();

    // only use first 10% of entries (for debugging)
    // uncomment to use
    totalEntries /= SCALEDOWN;

    for (int i=0; i<totalEntries; i++)
    {
        inputTree->GetEntry(i);

        if (rawEvent.chNo==0)
        {
            // treat DPP and waveform mode events differently
            if (rawEvent.evtType==1)
            {
                // DPP mode
                if (rawEvent.lgQ==65535)
                {
                    // ignore target changer events having an off-scale integrated
                    // charge - these are suspected to be retriggers and are NOT
                    // time-correlated with the facility RF time reference.
                    continue;
                }

                // assign the macropulse start time
                tcEvent.macroTime = (double)rawEvent.extTime*pow(2,32)+rawEvent.timetag;

                // Check to see if this is the start of a new DPP period
                // (i.e., a switch between DPP and waveform modes)
                if(prevEvtType==2)
                {
                    // the last event from this channel was in waveform mode
                    // so we must be in a new DPP mode period.
                    // Set evtNo = 1 to indicate this.
                    // When we sort through the other channel specific
                    // trees we'll know where the mode switches are. 
                    tcEvent.modeChange = 1;
                }

                else
                {
                    // NOT the start of a new DPP period (99.8% of the time)
                    tcEvent.modeChange = 0;
                }

                // Check for digitizer error (incrementing extTime before clearing timetag)
                if (rawEvent.extTime > extTimePrev && rawEvent.timetag > pow(2,32)-1000)
                {
                    error << "Found a target changer event with a timestamp-reset failure (i.e., extTime incremented before timetag was reset to 0). MacroNo = " << tcEvent.macroNo << ", extTime = " << rawEvent.extTime << ", extTimePrev = " << extTimePrev << ", timetag = " << rawEvent.timetag << ", timetagPrev = " << timetagPrev << endl;
                    error << "Skipping to next target changer event..." << endl;
                    continue;
                }

                tcEvent.targetPos = assignTargetPos(rawEvent.lgQ);

                /*************************************************************/
                // all the variables are updated: fill tree
                targetChangerTree->Fill();
                tcEvent.macroNo++;

                // before we pull the next event, save the previous event's data so we can
                // compare them with the new event's data
                extTimePrev = rawEvent.extTime;
                timetagPrev = rawEvent.timetag;

                if(tcEvent.macroNo%100==0)
                {
                    cout << "Target position = " << tcEvent.targetPos << ", macroNo = " << tcEvent.macroNo << "\r";
                    fflush(stdout);
                }
            }
            prevEvtType = rawEvent.evtType;
            // move to next event in the loop
        }
    }

    // loop finished
    // all target changer events sorted into targetChanger and targetChangerTreeW
    // time to move to detector channels

    cout << endl << "... done." << endl;

    // Now that macropulses have been assigned, point macropulse variables to
    // the target changer tree
    targetChangerTree->SetBranchAddress("macroNo",&tcEvent.macroNo);
    targetChangerTree->SetBranchAddress("macroTime",&tcEvent.macroTime);
    targetChangerTree->SetBranchAddress("modeChange",&tcEvent.modeChange);
    targetChangerTree->SetBranchAddress("targetPos",&tcEvent.targetPos);
}


// assign correct macropulse data (time reference of macropulse, target
// position) to each detector DPP event
void processDPPEvents(vector<TTree*> orchardRaw, vector<TTree*> orchardProcessed, TTree* targetChangerTree, ofstream& error)
{
    // skip channel 6 (scavenger) events
    for(int detIndex=2; detIndex<=4; detIndex+=2)
    {
        error << "Starting DPP processing on channel " << detIndex << endl;

        // Attach variables to correct channel-specific tree
        orchardRaw[detIndex/2]->SetBranchAddress("timetag",&rawEvent.timetag);
        orchardRaw[detIndex/2]->SetBranchAddress("extTime",&rawEvent.extTime);
        orchardRaw[detIndex/2]->SetBranchAddress("fineTime",&rawEvent.fineTime);
        orchardRaw[detIndex/2]->SetBranchAddress("sgQ",&rawEvent.sgQ);
        orchardRaw[detIndex/2]->SetBranchAddress("lgQ",&rawEvent.lgQ);
        orchardRaw[detIndex/2]->SetBranchAddress("waveform",&rawEvent.waveform);

        // reset event counter in preparation for looping through events
        procEvent.evtNo = 0; 

        // point to the first macropulse in preparation for assigning events to
        // macropulses
        targetChangerTree->GetEntry(0);

        // zero the dummy variables for comparing adjacent events
        double prevMacroTime = 0;
        double prevCompleteTime = 0;
        double extTimePrev = 0;
        double timetagPrev = 0;
        int prevEvtType = 0;

        // The experiment setup has a cable and electronics delay (TIME_OFFSET)
        // that is unique to each channel. All times are taken relative to the
        // macropulse state time.
        double TIME_OFFSET;

        switch(detIndex)
        {
            //case 0:
            //    TIME_OFFSET = 0;
            case 2:
                TIME_OFFSET = MACROPULSE_OFFSET;
                break;
            case 4:
                TIME_OFFSET = MACROPULSE_OFFSET;
                break;
            case 6:
                TIME_OFFSET = MACROPULSE_OFFSET+SCAVENGER_OFFSET;
                break;
            default:
                cout << "Error: non-existent channel index given by detIndex" << endl;
                exit(1);
        }

        // start looping through detector events and target changer events,
        // assigning detector events to macropulses.
        cout << endl << "Processing channel " << detIndex << " events..." << endl;

        long DPPTreeEntries = orchardRaw[detIndex/2]->GetEntries();
        long targetChangerTreeEntries = targetChangerTree->GetEntries();

        // Uncomment to take only first 10% of events (testing mode)
        DPPTreeEntries /= SCALEDOWN;

        for (int index=0; index<DPPTreeEntries; index++)
        {
            orchardRaw[detIndex/2]->GetEntry(index);

            // calculate this event's time
            procEvent.completeTime = (double)rawEvent.extTime*pow(2,32)+rawEvent.timetag+TIME_OFFSET;

            // Fine times are unused on channel 2 (monitor)
            if(detIndex!=2)
            {
                procEvent.completeTime+=rawEvent.fineTime*((double)2/1024);
            }

            // Ignore events that have problems with their time structure:
            // check to see if there was a timestamp-reset failure (the
            // extended timestamp incremented, but the 32-bit timestamp
            // didn't reset). If so, discard the event and move on.

            if (rawEvent.extTime > extTimePrev && rawEvent.timetag > pow(2,32)-10000)
            {
                error << endl << "MACRO NO " << tcEvent.macroNo << endl;
                error << "Found a detector event with a timestamp-reset failure (i.e., extTime incremented before timetag was reset to 0). completeTime = " << procEvent.completeTime << endl;
                error << "Skipping to next event..." << endl;
                continue;
            }

            /*****************************************************************/
            // Match events to their correct macropulses
            /*****************************************************************/

            // Next, check to see if this is event is the first one
            // of a new DPP mode period (with corresponding time reset).

            if (procEvent.completeTime < prevCompleteTime)
            {
                // new DPP mode period, triggered by this event's
                // timestamp wrapping back around to time zero

                // move to the next macropulse of the target changer
                error << "DPP/waveform mode switch at macroNo " << tcEvent.macroNo << endl;
                do
                {
                    if (tcEvent.macroNo+1>=targetChangerTreeEntries)
                    {
                        // ran out of new macropulses - start a series
                        // of breaks to exit the main loop
                        cout << endl << "Exceeded the last macropulse on targetChangerTree. Exiting loop." << endl;
                        break;
                    }

                    prevMacroTime = tcEvent.macroTime;
                    targetChangerTree->GetEntry(tcEvent.macroNo+1);
                    if(tcEvent.macroNo == 67387)
                    {
                        error << "macroNo = 67387, macroTime = " << tcEvent.macroTime << ", prevMacroTime = " << prevMacroTime << ", completeTime = " << procEvent.completeTime << ", prevCompleteTime = " << prevCompleteTime << endl;
                    }
                }
                while (prevMacroTime < tcEvent.macroTime);

                if (tcEvent.macroNo+1>=targetChangerTreeEntries)
                {
                    // exit the main loop
                    break;
                }

                // targetChangerTree is now pointing at the first
                // macropulse of the next DPP mode period
                //cout << "Start new DPP period @ macroNo " << tcEvent.macroNo << "\r";
                //fflush(stdout);

                // reset the event counter for this channel
                procEvent.evtNo = 0;
            }

            // Check to see if enough time has elapsed that the current
            // macropulse has expired (i.e., beam has gone off at end of macro).
            else if (procEvent.completeTime-tcEvent.macroTime-TIME_OFFSET > 0.95*MACRO_PERIOD)
            {
                // macropulse expired
                // time to switch to the next macropulse

                prevMacroTime = tcEvent.macroTime;

                // test to see if we've reached the last macropulse
                if (tcEvent.macroNo+1>=targetChangerTreeEntries)
                {
                    cout << endl << "Exceeded the last macropulse on targetChangerTree. Exiting loop." << endl;
                    break;
                }

                targetChangerTree->GetEntry(tcEvent.macroNo+1);

                //if(tcEvent.macroNo == 67387)
                //{
                //    error << "macroNo = 67387, macroTime = " << tcEvent.macroTime << ", prevMacroTime = " << prevMacroTime << ", completeTime = " << procEvent.completeTime << ", prevCompleteTime = " << prevCompleteTime << ", tc mode change = " << tcEvent.modeChange << endl;
                //}

                if (tcEvent.modeChange==1)
                {
                    // this macropulse is the first macropulse of a new DPP mode

                    // Move detector event index forward until it reaches a new DPP mode period
                    do
                    {
                        if (index==DPPTreeEntries)
                        {
                            cout << "Finished " << detIndex << " tree." << endl;
                            break;
                        }

                        index++;
                        orchardRaw[detIndex/2]->GetEntry(index);
                        prevCompleteTime = procEvent.completeTime;
                        procEvent.completeTime = (double)rawEvent.extTime*pow(2,32)+rawEvent.timetag+TIME_OFFSET;

                        if(detIndex!=2)
                        {
                            procEvent.completeTime+=rawEvent.fineTime*((double)2/1024);
                        }
                    }

                    while (prevCompleteTime < procEvent.completeTime/* || extTime > 0*/);

                    if (index==DPPTreeEntries)
                    {
                        break;
                    }

                    //error << "Start new DPP period @ macroNo " << tcEvent.macroNo << ", macroTime = " << tcEvent.macroTime << ", completeTime = " << procEvent.completeTime << endl;
                    //fflush(stdout);
                }

                // Lastly, check to see if there are any beam
                // anomalies (that is, the new macropulse is
                // out-of-sync with respect to the RF frequency )
                else if (!(tcEvent.macroTime-prevMacroTime > MACRO_PERIOD*(1-SYNC_WINDOW) &&
                           tcEvent.macroTime-prevMacroTime < MACRO_PERIOD*(1+SYNC_WINDOW)) &&

                         !(tcEvent.macroTime-prevMacroTime > 2*(MACRO_PERIOD*(1-SYNC_WINDOW)) &&
                           tcEvent.macroTime-prevMacroTime < 2*(MACRO_PERIOD*(1+SYNC_WINDOW))) &&

                         !(tcEvent.macroTime-prevMacroTime > 3*(MACRO_PERIOD*(1-SYNC_WINDOW)) &&
                           tcEvent.macroTime-prevMacroTime < 3*(MACRO_PERIOD*(1+SYNC_WINDOW))) &&

                         !(tcEvent.macroTime-prevMacroTime > 4*(MACRO_PERIOD*(1-SYNC_WINDOW)) &&
                           tcEvent.macroTime-prevMacroTime < 4*(MACRO_PERIOD*(1+SYNC_WINDOW)))
                        )
                {
                    // there is a beam anomaly - cycle through
                    // macropulses until we find another macropulse
                    // in sync w/ RF period

                    double beamAnomalyStart = tcEvent.macroTime;

                    error << endl << "MACRO NO " << tcEvent.macroNo << endl;
                    error << "Macropulse timing anomaly found at macroTime = " << beamAnomalyStart << endl;
                    error << "macroTime - prevMacroTime = " << tcEvent.macroTime-prevMacroTime << endl;
                    error << "Looping through macropulse list..." << endl;

                    do
                    {
                        if (tcEvent.macroNo+1 >= targetChangerTreeEntries)
                        {
                            cout << endl << "Exceeded the last macropulse on targetChangerTree. Exiting loop." << endl;
                            break;
                        }

                        prevMacroTime = tcEvent.macroTime;
                        targetChangerTree->GetEntry(tcEvent.macroNo+1);
                        error << "macroTime = " << tcEvent.macroTime << "; macroTime-prevMacroTime = " << tcEvent.macroTime-prevMacroTime << endl;
                    } while (!(tcEvent.macroTime-prevMacroTime > MACRO_PERIOD*(1-SYNC_WINDOW) &&
                               tcEvent.macroTime-prevMacroTime < MACRO_PERIOD*(1+SYNC_WINDOW)) &&

                             !(tcEvent.macroTime-prevMacroTime > 2*(MACRO_PERIOD*(1-SYNC_WINDOW)) &&
                               tcEvent.macroTime-prevMacroTime < 2*(MACRO_PERIOD*(1+SYNC_WINDOW))) &&

                             !(tcEvent.macroTime-prevMacroTime > 3*(MACRO_PERIOD*(1-SYNC_WINDOW)) &&
                               tcEvent.macroTime-prevMacroTime < 3*(MACRO_PERIOD*(1+SYNC_WINDOW))) &&

                             !(tcEvent.macroTime-prevMacroTime > 4*(MACRO_PERIOD*(1-SYNC_WINDOW)) &&
                               tcEvent.macroTime-prevMacroTime < 4*(MACRO_PERIOD*(1+SYNC_WINDOW)))
                            );

                    error << "Beam anomaly finished; new macroTime is " << tcEvent.macroTime << ", macroNo " << tcEvent.macroNo << ", completeTime = " << procEvent.completeTime << endl;

                    // We've found the next 'good' macropulse.
                    // Need to move forward in the raw tree to point at the
                    // same area of good beam

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

                    if (beamAnomalyStart > tcEvent.macroTime)
                    {
                        // beam anomaly DID span a DPP/waveform mode
                        // move the raw tree forward past mode change
                        error << "Beam anomaly spanned DPP/waveform mode change; moving the channel index forward to the next DPP mode period." << endl;
                        do
                        {
                            if (index==DPPTreeEntries)
                            {
                                error << "Reached the end of the raw tree; end looping." << endl;
                                break;
                            }

                            index++;
                            orchardRaw[detIndex/2]->GetEntry(index);
                            prevCompleteTime = procEvent.completeTime;
                            procEvent.completeTime = (double)rawEvent.extTime*pow(2,32)+rawEvent.timetag+TIME_OFFSET;
                        }

                        // keep looping until the timestamps reset to 0 for
                        // this channel (i.e., new DPP period)
                        while (prevCompleteTime < procEvent.completeTime);

                        // if we've reached the end of the raw tree events,
                        // exit the loop
                        if (index==DPPTreeEntries)
                        {
                            break;
                        }

                        error << "Finished moving channel index forward to the next DPP mode period." << endl;
                    }

                    // continue moving index forward until we reach the
                    // current macropulse time
                    do
                    {
                        index++;
                        orchardRaw[detIndex/2]->GetEntry(index);
                        procEvent.completeTime = (double)rawEvent.extTime*pow(2,32)+rawEvent.timetag+TIME_OFFSET;
                    }
                    while (procEvent.completeTime < tcEvent.macroTime);
                }

                // Lastly, because we're in a new macropulse: reset the event counter
                procEvent.evtNo = 0;
            }

            /*****************************************************************/
            // Add processed events into destination tree
            /*****************************************************************/

            // Assign the macropulse variables for the event
            procEvent.macroNo = tcEvent.macroNo;
            procEvent.macroTime = tcEvent.macroTime;
            procEvent.targetPos = tcEvent.targetPos;
            procEvent.sgQ = rawEvent.sgQ;
            procEvent.lgQ = rawEvent.lgQ;
            //procEvent.waveform = rawEvent.waveform;

            // only add events that come while beam is on, during the macropulse 
            double timeDiff = procEvent.completeTime-tcEvent.macroTime;

            if (timeDiff > 0 && timeDiff < MACROPULSE_LENGTH)
            {
                fillProcessedTree(orchardProcessed[detIndex/2]);
            }

            // in preparation for looping to the next event, update counters
            numberOfDPPs++;
            procEvent.evtNo++;

            prevCompleteTime = procEvent.completeTime;
            prevEvtType = rawEvent.evtType;
            extTimePrev = rawEvent.extTime;
            timetagPrev = rawEvent.timetag;

            if(index%10000==0)
            {
                cout << "Event " << index << "\r";
                fflush(stdout);
            }
        }
    }
}

// assign target position to each waveform event
void processWaveformEvents(vector<TTree*> orchardRawW, vector<TTree*> orchardProcessedW, TTree* targetChangerTree, ofstream& error)
{
    for(int detIndex=0; detIndex<=4; detIndex+=2)
    {
        // Attach variables to correct waveform tree
        orchardRawW[detIndex/2]->SetBranchAddress("timetag",&rawEvent.timetag);
        orchardRawW[detIndex/2]->SetBranchAddress("extTime",&rawEvent.extTime);
        orchardRawW[detIndex/2]->SetBranchAddress("evtNo",&rawEvent.evtNo);
        orchardRawW[detIndex/2]->SetBranchAddress("waveform",&rawEvent.waveform);

        int totalEntries = orchardRawW[detIndex/2]->GetEntries();

        int prevMacroNo = 0;

        targetChangerTree->GetEntry(0);

        for(int i=0; i<totalEntries; i++)
        {
            orchardRawW[detIndex/2]->GetEntry(i);

            //error << "waveform evtNo = " << rawEvent.evtNo << endl;
            if(rawEvent.evtNo==0)
            {
                // Start of new waveform mode period:
                // Move forward to the next mode change in the target changer

                //error << "modeChange == " << tcEvent.modeChange << endl;
                do
                {
                    targetChangerTree->GetEntry(tcEvent.macroNo+1);
                //    error << "cycling through tcTree, macroNo = " << tcEvent.macroNo << endl;

                } while(tcEvent.modeChange==0 && tcEvent.macroNo+1<targetChangerTree->GetEntries());

                if(tcEvent.macroNo+1>=targetChangerTree->GetEntries())
                {
                    // reached end of targetChangerTree
                    break;
                }

                prevMacroNo=tcEvent.macroNo;
            }

            // calculate this event's time
            procEvent.completeTime = (double)rawEvent.extTime*pow(2,32)+rawEvent.timetag;

            // assign remaining event variables to processed event
            procEvent.macroNo = tcEvent.macroNo;
            procEvent.macroTime = tcEvent.macroTime;
            procEvent.targetPos = tcEvent.targetPos;

            error << "targ position = " << tcEvent.targetPos << endl;

            procEvent.evtNo = rawEvent.evtNo;
            procEvent.waveform = rawEvent.waveform;

            fillProcessedTree(orchardProcessedW[detIndex/2]);
            numberOfWaveforms++;
        }
    }
}

/******************************************************************************/

int main(int argc, char* argv[])
{
    /**************************************************************************
                        Check to see if a sort is necessary
    **************************************************************************/

    cout << endl << "Entering re-sort..." << endl;

    string runDir = argv[1];
    string runNo = argv[2];
    string outpath = argv[3];

    TFile *fileOut;

    stringstream treeName;
    treeName << "run" << runDir << "-" << runNo; 

    stringstream fileOutName;
    fileOutName << outpath <<"/analysis/run" << runDir << "/" << treeName.str() << "_sorted.root";
    fileOut = new TFile(fileOutName.str().c_str(),"UPDATE");

    // Check to see if sorted trees already exist in the output file. If they
    // do, exit - our job is already done.
    if(fileOut->Get("targetChangerTree"))
    {
        cout << fileOutName.str() << " already exists. Skipping ./resort..." << endl;
        exit(0);
    }

    /***************************************************************************
                                Prepare for sorting
    ***************************************************************************/

    // Create holders for the channel-specific trees, separated by event type and stage of
    // processing stage (i.e., whether events in the tree are unassigned to macropulses (raw) or
    // assigned to macropulses (processed))
    vector<TTree*> orchardRaw;       // channel-specific DPP events NOT assigned to macropulses
    vector<TTree*> orchardProcessed; // channel-specific DPP events assigned to macropulses
    vector<TTree*> orchardRawW;      // channel-specific waveform events NOT assigned to macropulses
    vector<TTree*> orchardProcessedW;// channel-specific waveform events assigned to macropulses

    // Create the new empty trees
    // Each channel has a separate tree for DPP data and for waveform mode data
    for(int i=0; (size_t)i<activeDPPChannels.size(); i++)
    {
        orchardRaw.push_back(new TTree((activeDPPChannels[i]+"RawTree").c_str(),""));
        branchRaw(orchardRaw[i]);

        orchardProcessed.push_back(new TTree((activeDPPChannels[i]+"ProcessedTree").c_str(),""));
        branchProc(orchardProcessed[i]);
    }

    for(int i=0; (size_t)i<activeWaveformChannels.size(); i++)
    {
        orchardRawW.push_back(new TTree((activeWaveformChannels[i]+"RawTreeW").c_str(),""));
        branchRawW(orchardRawW[i]);

        orchardProcessedW.push_back(new TTree((activeWaveformChannels[i]+"ProcessedTreeW").c_str(),""));
        branchProcW(orchardProcessedW[i]);
    }

    // Create tree for tracking macropulses in this sub-run
    TTree* targetChangerTree = new TTree("targetChangerTree","");
    branchTargetChanger(targetChangerTree);

    /***************************************************************************
                               Prepare to sort input file
    ***************************************************************************/

    // Set cout precision to display entire timestamp values
    cout.precision(13);

    // Access the input file containing the raw events tree
    TFile *fileIn;
    stringstream fileInName;
    fileInName << outpath <<"/analysis/run" << runDir << "/" << treeName.str() << "_raw.root";
    fileIn = new TFile(fileInName.str().c_str(),"READ");

    if(!fileIn->Get("tree"))
    {
        cout << "Error: failed to find raw tree in " << fileInName.str() << endl;
        exit(1);
    }

    // Access the input file tree
    TTree *inputTree;
    inputTree = (TTree*)fileIn->Get("tree");

    // Create an error log where sorting errors can be recorded for review
    ofstream error;
    stringstream errorName;
    errorName << outpath <<  "/analysis/run" << runDir << "/" << treeName.str() << "_error.log";
    error.open(errorName.str());
    error.precision(13);

    // link the tree from the input file to our event variables
    inputTree->SetBranchAddress("chNo",&rawEvent.chNo);
    inputTree->SetBranchAddress("evtType",&rawEvent.evtType);
    inputTree->SetBranchAddress("timetag",&rawEvent.timetag);
    inputTree->SetBranchAddress("extTime",&rawEvent.extTime);
    inputTree->SetBranchAddress("sgQ",&rawEvent.sgQ);
    inputTree->SetBranchAddress("lgQ",&rawEvent.lgQ);
    inputTree->SetBranchAddress("fineTime",&rawEvent.fineTime);
    inputTree->SetBranchAddress("waveform",&rawEvent.waveform);

    /***************************************************************************
                            Sort all events from input file
    ***************************************************************************/

    // Separate events
    // from the input tree into channel-specific trees
    separateByChannel(inputTree, orchardRaw, orchardRawW);

    // Next, extract target changer events from the input tree and add to the
    // target changer trees (DPP and waveform), assigning a macropulse to each
    // target changer event. 
    processTargetChanger(inputTree, targetChangerTree, error);

    // Last, now that the macropulse structure is assigned by the target changer
    // events, we can assign detector events to the correct macropulse.
    processDPPEvents(orchardRaw, orchardProcessed, targetChangerTree, error);
    processWaveformEvents(orchardRawW, orchardProcessedW, targetChangerTree, error);

    cout << "Total number of DPP-mode events processed = " << numberOfDPPs << endl;
    cout << "Total number of waveform-mode events processed = " << numberOfWaveforms << endl;
    /*cout << "Total number of ch0 waveform-mode events processed = " << numberOfCh0Waveforms << endl;
    cout << "Total number of ch2 waveform-mode events processed = " << numberOfCh2Waveforms << endl;
    cout << "Total number of ch4 waveform-mode events processed = " << numberOfCh4Waveforms << endl;
    */

    /***************************************************************************
                                 Clean up and exit
    ***************************************************************************/

    // Delete raw trees
    for(vector<TTree*>::iterator orchardRawIterator = orchardRaw.begin();
            orchardRawIterator != orchardRaw.end();
            orchardRawIterator++)
    {
        delete *orchardRawIterator;
    }

    for(vector<TTree*>::iterator orchardRawWIterator = orchardRawW.begin();
            orchardRawWIterator != orchardRawW.end();
            orchardRawWIterator++)
    {
        delete *orchardRawWIterator;
    }

    /*gDirectory->Delete("ch0RawTree;*");
    gDirectory->Delete("ch2RawTree;*");
    gDirectory->Delete("ch4RawTree;*");
    gDirectory->Delete("ch6RawTree;*");
    gDirectory->Delete("ch0RawTreeW;*");
    gDirectory->Delete("ch0RawTreeW;*");
    gDirectory->Delete("ch0RawTreeW;*");
    */

    fileOut->Write();
    fileOut->Close();
}
