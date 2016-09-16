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

#include "../include/dataStructures.h"
#include "../include/analysisConstants.h"
#include "../include/physicalConstants.h"
#include "../include/resort.h"
#include "../include/branches.h"

using namespace std;

extern RawEvent rawEvent;
extern SortedEvent sortedEvent;
extern ProcessedEvent procEvent;
extern TargetChangerEvent tcEvent;

long sortedNumberOfCh0Waveforms = 0;
long sortedNumberOfCh2Waveforms = 0;
long sortedNumberOfCh4Waveforms = 0;

long sortedNumberOfDPPs = 0;
long sortedNumberOfWaveforms = 0;

const double DEBUG_SCALEDOWN = 1; // (for debugging) only sort (total/DEBUG_SCALEDOWN) events

void fillRawTree(TTree* tree)
{
    // To save only 1 DPP event's waveform in every 10000 events, uncomment this
    // (Note - will reduce output tree size by ~half)
    /*if ((int)completeTime%10000 != 0 && evtType==1)
      {
      waveform->clear(); 
      }*/

    tree->Fill();

    if(sortedEvent.evtType==2)
    {
        if(sortedEvent.chNo==0)
        {
            sortedNumberOfCh0Waveforms++;
        }

        if(sortedEvent.chNo==2)
        {
            sortedNumberOfCh2Waveforms++;
        }

        if(sortedEvent.chNo==4)
        {
            sortedNumberOfCh4Waveforms++;
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
    sortedEvent.waveform->clear();

    // fill waveform with wavelet data from an event
    // (the raw tree's branch points to dummyWaveform)
    for(int k=0; (size_t)k<dummyWaveform->size(); k++)
    {
        sortedEvent.waveform->push_back(dummyWaveform->at(k));
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
void separateByChannel(string rawFileName, string tempFileName, vector<TTree*>& orchardRaw, vector<TTree*>& orchardRawW)
{
    /*sortedNumberOfCh0Waveforms = 0;
    sortedNumberOfCh2Waveforms = 0;
    sortedNumberOfCh4Waveforms = 0;

    sortedNumberOfDPPs = 0;
    sortedNumberOfWaveforms = 0;
    */

    TFile* tempFile = new TFile(tempFileName.c_str(), "RECREATE");

    // Create the new empty trees
    // Each channel has a separate tree for DPP data and for waveform mode data
    for(int i=0; (size_t)i<activeDPPChannels.size(); i++)
    {
        orchardRaw.push_back(new TTree((activeDPPChannels[i]+"RawTree").c_str(),""));
        branchRaw(orchardRaw[i]);
        orchardRaw[i]->SetDirectory(tempFile);
    }

    for(int i=0; (size_t)i<activeWaveformChannels.size(); i++)
    {
        orchardRawW.push_back(new TTree((activeWaveformChannels[i]+"RawTreeW").c_str(),""));
        branchRawW(orchardRawW[i]);
        orchardRawW[i]->SetDirectory(tempFile);
    }

    cout << "Separating events by channel and event type..." << endl;

    TFile* rawFile = new TFile(rawFileName.c_str(),"READ");
    TTree* inputTree = (TTree*)rawFile->Get("tree");

    if(!rawFile->Get("tree"))
    {
        cout << "Error: failed to find raw tree in " << rawFileName << endl;
        exit(1);
    }

    // link the tree from the input file to our event variables
    inputTree->SetBranchAddress("chNo",&sortedEvent.chNo);
    inputTree->SetBranchAddress("evtType",&sortedEvent.evtType);
    inputTree->SetBranchAddress("timetag",&sortedEvent.timetag);
    inputTree->SetBranchAddress("extTime",&sortedEvent.extTime);
    inputTree->SetBranchAddress("sgQ",&sortedEvent.sgQ);
    inputTree->SetBranchAddress("lgQ",&sortedEvent.lgQ);
    inputTree->SetBranchAddress("fineTime",&sortedEvent.fineTime);
    inputTree->SetBranchAddress("waveform",&sortedEvent.waveform);

    int prevEvtType = 0;
    
    int totalEntries = inputTree->GetEntries();

    totalEntries /= DEBUG_SCALEDOWN; // for debugging:
                               // use to sort only a subset of total events

    // To uniquely identify each event, we assign each channel's events an event
    // number (evtNo), which is the event's order in its macropulse.
    vector<int> evtNo(activeDPPChannels.size());
    
    // Loop through all events in input tree and separate them into channel-
    // specific trees
    for (int index=0; index<totalEntries; index++)
    {
        inputTree->GetEntry(index);

        // if the event is a target changer event (chNo = 0),
        // processTargetChanger will handle it separately because target
        // changer events are special, non-detector events
        /*if(sortedEvent.chNo==0)
        {
            continue;
        }*/

        if (sortedEvent.evtType==1)
        {
            // DPP mode

            fillRawTree(orchardRaw[sortedEvent.chNo/2]);
        }

        if (sortedEvent.evtType==2)
        {
            // Waveform mode
            // Because waveform mode events are not time referenced to target
            // changer events, we need to manually increment the evtNo here 
            if (prevEvtType==1)
            {
                // New waveform mode period
                for(auto &value: evtNo)
                {
                    value = 0;
                }
            }

            sortedEvent.evtNo = evtNo[sortedEvent.chNo/2];
            fillRawTree(orchardRawW[sortedEvent.chNo/2]);
            evtNo[sortedEvent.chNo/2]++;
        }

        if(index%10000==0)
        {
            cout << "processed " << index << " events.\r";
            fflush(stdout);
        }

        prevEvtType = sortedEvent.evtType;
    }

    cout << "Separated " << totalEntries << " events into channels 2, 4 and 6." << endl;
    tempFile->Write();

    rawFile->Close();
}

// assign a macropulse, target position, etc to every target changer event
void processTargetChanger(string rawFileName, TFile*& sortedFile, ofstream& error)
{
    TTree* targetChangerTree = new TTree("targetChangerTree","");
    branchTargetChanger(targetChangerTree);
    targetChangerTree->SetDirectory(sortedFile);

    TFile* rawData = new TFile(rawFileName.c_str(),"READ");
    TTree* inputTree = (TTree*)rawData->Get("tree");

    inputTree->SetBranchAddress("chNo",&sortedEvent.chNo);
    inputTree->SetBranchAddress("evtType",&sortedEvent.evtType);
    inputTree->SetBranchAddress("timetag",&sortedEvent.timetag);
    inputTree->SetBranchAddress("extTime",&sortedEvent.extTime);
    inputTree->SetBranchAddress("sgQ",&sortedEvent.sgQ);
    inputTree->SetBranchAddress("lgQ",&sortedEvent.lgQ);
    inputTree->SetBranchAddress("fineTime",&sortedEvent.fineTime);
    inputTree->SetBranchAddress("waveform",&sortedEvent.waveform);

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
    totalEntries /= DEBUG_SCALEDOWN;

    for(int i=0; i<totalEntries; i++)
    {
        inputTree->GetEntry(i);

        if(sortedEvent.chNo==0)
        {
            // treat DPP and waveform mode events differently
            if(sortedEvent.evtType==1)
            {
                // DPP mode
                if (sortedEvent.lgQ==65535)
                {
                    // ignore target changer events having an off-scale integrated
                    // charge - these are suspected to be retriggers and are NOT
                    // time-correlated with the facility RF time reference.
                    continue;
                }

                // assign the macropulse start time
                tcEvent.macroTime = (double)sortedEvent.extTime*pow(2,32)+sortedEvent.timetag;

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
                if (sortedEvent.extTime > extTimePrev && sortedEvent.timetag > pow(2,32)-1000)
                {
                    error << "Found a target changer event with a timestamp-reset failure (i.e., extTime incremented before timetag was reset to 0). MacroNo = " << tcEvent.macroNo << ", extTime = " << sortedEvent.extTime << ", extTimePrev = " << extTimePrev << ", timetag = " << sortedEvent.timetag << ", timetagPrev = " << timetagPrev << endl;
                    error << "Skipping to next target changer event..." << endl;
                    continue;
                }

                tcEvent.targetPos = assignTargetPos(sortedEvent.lgQ);

                /*************************************************************/
                // all the variables are updated: fill tree
                targetChangerTree->Fill();
                tcEvent.macroNo++;

                // before we pull the next event, save the previous event's data so we can
                // compare them with the new event's data
                extTimePrev = sortedEvent.extTime;
                timetagPrev = sortedEvent.timetag;

                if(tcEvent.macroNo%100==0)
                {
                    cout << "Target position = " << tcEvent.targetPos << ", macroNo = " << tcEvent.macroNo << "\r";
                    fflush(stdout);
                }
            }
            prevEvtType = sortedEvent.evtType;
            // move to next event in the loop
        }
    }

    // loop finished
    // all target changer events sorted into targetChanger and targetChangerTreeW
    // time to move to detector channels

    cout << endl << "... done." << endl;

    sortedFile->cd();
    targetChangerTree->Write();
}


// assign correct macropulse data (time reference of macropulse, target
// position) to each detector DPP event
void processDPPEvents(TFile*& sortedFile, vector<TTree*>& orchardRaw, vector<TTree*>& orchardProcessed, ofstream& error)
{
    TTree* targetChangerTree = (TTree*)sortedFile->Get("targetChangerTree");

    // Now that macropulses have been assigned, point macropulse variables to
    // the target changer tree
    targetChangerTree->SetBranchAddress("macroNo",&tcEvent.macroNo);
    targetChangerTree->SetBranchAddress("macroTime",&tcEvent.macroTime);
    targetChangerTree->SetBranchAddress("modeChange",&tcEvent.modeChange);
    targetChangerTree->SetBranchAddress("targetPos",&tcEvent.targetPos);

    // Create the new empty trees
    // Each channel has a separate tree for DPP data and for waveform mode data
    for(int i=0; (size_t)i<activeDPPChannels.size(); i++)
    {
        orchardProcessed.push_back(new TTree((activeDPPChannels[i]+"ProcessedTree").c_str(),""));
        branchProc(orchardProcessed[i]);
        orchardProcessed[i]->SetDirectory(sortedFile);
    }

    // skip channel 6 (scavenger) events
    for(int detIndex=2; detIndex<=4; detIndex+=2)
    {
        error << "Starting DPP processing on channel " << detIndex << endl;

        // Attach variables to correct channel-specific tree
        orchardRaw[detIndex/2]->SetBranchAddress("timetag",&sortedEvent.timetag);
        orchardRaw[detIndex/2]->SetBranchAddress("extTime",&sortedEvent.extTime);
        orchardRaw[detIndex/2]->SetBranchAddress("fineTime",&sortedEvent.fineTime);
        orchardRaw[detIndex/2]->SetBranchAddress("sgQ",&sortedEvent.sgQ);
        orchardRaw[detIndex/2]->SetBranchAddress("lgQ",&sortedEvent.lgQ);
        orchardRaw[detIndex/2]->SetBranchAddress("waveform",&sortedEvent.waveform);

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

        for (int index=0; index<DPPTreeEntries; index++)
        {
            orchardRaw[detIndex/2]->GetEntry(index);

            // calculate this event's time
            procEvent.completeTime = (double)sortedEvent.extTime*pow(2,32)+sortedEvent.timetag+TIME_OFFSET;

            // Fine times are unused on channel 2 (monitor)
            if(detIndex!=2)
            {
                procEvent.completeTime+=sortedEvent.fineTime*((double)2/1024);
            }

            // Ignore events that have problems with their time structure:
            // check to see if there was a timestamp-reset failure (the
            // extended timestamp incremented, but the 32-bit timestamp
            // didn't reset). If so, discard the event and move on.

            if (sortedEvent.extTime > extTimePrev && sortedEvent.timetag > pow(2,32)-10000)
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
                    if(tcEvent.macroNo == 1320)
                    {
                        error << "macroNo = 1320, macroTime = " << tcEvent.macroTime << ", prevMacroTime = " << prevMacroTime << ", completeTime = " << procEvent.completeTime << ", prevCompleteTime = " << prevCompleteTime << endl;
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
            else if (procEvent.completeTime-tcEvent.macroTime-TIME_OFFSET > MACRO_LENGTH)
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
                        procEvent.completeTime = (double)sortedEvent.extTime*pow(2,32)+sortedEvent.timetag+TIME_OFFSET;

                        if(detIndex!=2)
                        {
                            procEvent.completeTime+=sortedEvent.fineTime*((double)2/1024);
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
                            procEvent.completeTime = (double)sortedEvent.extTime*pow(2,32)+sortedEvent.timetag+TIME_OFFSET;
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
                        procEvent.completeTime = (double)sortedEvent.extTime*pow(2,32)+sortedEvent.timetag+TIME_OFFSET;
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
            procEvent.sgQ = sortedEvent.sgQ;
            procEvent.lgQ = sortedEvent.lgQ;
            procEvent.waveform = sortedEvent.waveform;

            // only add events that come while beam is on, during the macropulse 
            double timeDiff = procEvent.completeTime-tcEvent.macroTime;

            if (timeDiff > 0 && timeDiff < MACRO_LENGTH)
            {
                fillProcessedTree(orchardProcessed[detIndex/2]);
            }

            // in preparation for looping to the next event, update counters
            sortedNumberOfDPPs++;
            procEvent.evtNo++;

            prevCompleteTime = procEvent.completeTime;
            prevEvtType = sortedEvent.evtType;
            extTimePrev = sortedEvent.extTime;
            timetagPrev = sortedEvent.timetag;

            if(index%10000==0)
            {
                cout << "Event " << index << "\r";
                fflush(stdout);
            }
        }
    }
    cout << "Total number of DPP-mode events processed = " << sortedNumberOfDPPs << endl;
}

// assign target position to each waveform event
void processWaveformEvents(TFile*& sortedFile, vector<TTree*>& orchardRawW, vector<TTree*>& orchardProcessedW, ofstream& error)
{
    TTree* targetChangerTree = (TTree*)sortedFile->Get("targetChangerTree");

    // Now that macropulses have been assigned, point macropulse variables to
    // the target changer tree
    targetChangerTree->SetBranchAddress("macroNo",&tcEvent.macroNo);
    targetChangerTree->SetBranchAddress("macroTime",&tcEvent.macroTime);
    targetChangerTree->SetBranchAddress("modeChange",&tcEvent.modeChange);
    targetChangerTree->SetBranchAddress("targetPos",&tcEvent.targetPos);



    for(int i=0; (size_t)i<activeWaveformChannels.size(); i++)
    {
        orchardProcessedW.push_back(new TTree((activeWaveformChannels[i]+"ProcessedTreeW").c_str(),""));
        branchProcW(orchardProcessedW[i]);
        //orchardProcessedW[i]->SetDirectory(tempFile);
    }

    for(int detIndex=0; detIndex<=4; detIndex+=2)
    {
        // Attach variables to correct waveform tree
        orchardRawW[detIndex/2]->SetBranchAddress("timetag",&sortedEvent.timetag);
        orchardRawW[detIndex/2]->SetBranchAddress("extTime",&sortedEvent.extTime);
        orchardRawW[detIndex/2]->SetBranchAddress("evtNo",&sortedEvent.evtNo);
        orchardRawW[detIndex/2]->SetBranchAddress("waveform",&sortedEvent.waveform);

        int totalEntries = orchardRawW[detIndex/2]->GetEntries();

        int prevMacroNo = 0;

        targetChangerTree->GetEntry(0);

        for(int i=0; i<totalEntries; i++)
        {
            orchardRawW[detIndex/2]->GetEntry(i);

            //error << "waveform evtNo = " << sortedEvent.evtNo << endl;
            if(sortedEvent.evtNo==0)
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
            procEvent.completeTime = (double)sortedEvent.extTime*pow(2,32)+sortedEvent.timetag;

            // assign remaining event variables to processed event
            procEvent.macroNo = tcEvent.macroNo;
            procEvent.macroTime = tcEvent.macroTime;
            procEvent.targetPos = tcEvent.targetPos;

            error << "targ position = " << tcEvent.targetPos << endl;

            procEvent.evtNo = sortedEvent.evtNo;
            procEvent.waveform = sortedEvent.waveform;

            fillProcessedTree(orchardProcessedW[detIndex/2]);
            sortedNumberOfWaveforms++;
        }
    }
    cout << "Total number of waveform-mode events processed = " << sortedNumberOfWaveforms << endl;
}