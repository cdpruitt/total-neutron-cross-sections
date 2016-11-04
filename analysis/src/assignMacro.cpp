/******************************************************************************
                               assignMacro.cpp 
******************************************************************************/
// assignMacro.cpp reads detector and macropulse events from ROOT trees
// (produced by ./separate) and assigns macropulse information to each detector event.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"

#include "../include/dataStructures.h"
#include "../include/analysisConstants.h"
#include "../include/physicalConstants.h"
#include "../include/separate.h"
#include "../include/branches.h"

extern SeparatedEvent separatedEvent;
extern ProcessedEvent procEvent;
extern TargetChangerEvent tcEvent;

using namespace std;

void processDPPEvents(TFile*& sortedFile, vector<TTree*>& orchardRaw, vector<TTree*>& orchardProcessed)
{
    TTree* targetChangerTree = (TTree*)sortedFile->Get("targetChangerTree");

    // Now that macropulses have been assigned, point macropulse variables to
    // the target changer tree
    setBranchesProcessedTC(targetChangerTree);
        // Create the new empty trees
    // Each channel has a separate tree for DPP data and for waveform mode data
    for(int i=0; (size_t)i<activeDPPChannels.size(); i++)
    {
        orchardProcessed.push_back(new TTree((activeDPPChannels[i]+"ProcessedTree").c_str(),""));
        branchProc(orchardProcessed[i]);
        orchardProcessed[i]->SetDirectory(sortedFile);
    }

    long separatedNumberOfDPPs = 0;

    for(int detIndex=2; detIndex<NUMBER_OF_CHANNELS*2; detIndex+=2)
    {
        cerr << "Starting DPP processing on channel " << detIndex << endl;

        // Attach variables to correct channel-specific tree
        setBranchesSeparated(orchardRaw[detIndex/2]);
        
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
                TIME_OFFSET = MACROPULSE_OFFSET-VETO_OFFSET;
                break;
            default:
                cerr << "Error: non-existent channel index given by detIndex" << endl;
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
            procEvent.completeTime = (double)separatedEvent.extTime*pow(2,32)+separatedEvent.timetag+TIME_OFFSET;

            // Fine times are unused on channel 2 (monitor)
            if(detIndex!=2)
            {
                procEvent.completeTime+=separatedEvent.fineTime*((double)2/1024);
            }

            // Ignore events that have problems with their time structure:
            // check to see if there was a timestamp-reset failure (the
            // extended timestamp incremented, but the 32-bit timestamp
            // didn't reset). If so, discard the event and move on.

            if (separatedEvent.extTime > extTimePrev && separatedEvent.timetag > pow(2,32)-10000)
            {
                cerr << endl << "MACRO NO " << tcEvent.macroNo << endl;
                cerr << "Found a detector event with a timestamp-reset failure (i.e., extTime incremented before timetag was reset to 0). completeTime = " << procEvent.completeTime << endl;
                cerr << "Skipping to next event..." << endl;
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
                cout << "DPP/waveform mode switch at macroNo " << tcEvent.macroNo << "\r";
                fflush(stdout);

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

                if (tcEvent.modeChange==1)
                {
                    // this macropulse is the first macropulse of a new DPP mode

                    cout << "modeChange = 1; new DPP mode at macroNo " << tcEvent.macroNo << "\r";
                    fflush(stdout);

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
                        procEvent.completeTime = (double)separatedEvent.extTime*pow(2,32)+separatedEvent.timetag+TIME_OFFSET;

                        if(detIndex!=2)
                        {
                            procEvent.completeTime+=separatedEvent.fineTime*((double)2/1024);
                        }
                    }

                    while (prevCompleteTime < procEvent.completeTime/* || extTime > 0*/);

                    if (index==DPPTreeEntries)
                    {
                        break;
                    }
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

                    cerr << endl << "MACRO NO " << tcEvent.macroNo << endl;
                    cerr << "Macropulse timing anomaly found at macroTime = " << beamAnomalyStart << endl;
                    cerr << "macroTime - prevMacroTime = " << tcEvent.macroTime-prevMacroTime << endl;
                    cerr << "Looping through macropulse list..." << endl;

                    do
                    {
                        if (tcEvent.macroNo+1 >= targetChangerTreeEntries)
                        {
                            cout << endl << "Exceeded the last macropulse on targetChangerTree. Exiting loop." << endl;
                            break;
                        }

                        prevMacroTime = tcEvent.macroTime;
                        targetChangerTree->GetEntry(tcEvent.macroNo+1);
                        cerr << "macroTime = " << tcEvent.macroTime << "; macroTime-prevMacroTime = " << tcEvent.macroTime-prevMacroTime << endl;
                    } while (!(tcEvent.macroTime-prevMacroTime > MACRO_PERIOD*(1-SYNC_WINDOW) &&
                               tcEvent.macroTime-prevMacroTime < MACRO_PERIOD*(1+SYNC_WINDOW)) &&

                             !(tcEvent.macroTime-prevMacroTime > 2*(MACRO_PERIOD*(1-SYNC_WINDOW)) &&
                               tcEvent.macroTime-prevMacroTime < 2*(MACRO_PERIOD*(1+SYNC_WINDOW))) &&

                             !(tcEvent.macroTime-prevMacroTime > 3*(MACRO_PERIOD*(1-SYNC_WINDOW)) &&
                               tcEvent.macroTime-prevMacroTime < 3*(MACRO_PERIOD*(1+SYNC_WINDOW))) &&

                             !(tcEvent.macroTime-prevMacroTime > 4*(MACRO_PERIOD*(1-SYNC_WINDOW)) &&
                               tcEvent.macroTime-prevMacroTime < 4*(MACRO_PERIOD*(1+SYNC_WINDOW)))
                            );

                    cerr << "Beam anomaly finished; new macroTime is " << tcEvent.macroTime << ", macroNo " << tcEvent.macroNo << ", completeTime = " << procEvent.completeTime << endl;

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
                        cerr << "Beam anomaly spanned DPP/waveform mode change; moving the channel index forward to the next DPP mode period." << endl;
                        do
                        {
                            if (index==DPPTreeEntries)
                            {
                                cerr << "Reached the end of the raw tree; end looping." << endl;
                                break;
                            }

                            index++;
                            orchardRaw[detIndex/2]->GetEntry(index);
                            prevCompleteTime = procEvent.completeTime;
                            procEvent.completeTime = (double)separatedEvent.extTime*pow(2,32)+separatedEvent.timetag+TIME_OFFSET;
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

                        cerr << "Finished moving channel index forward to the next DPP mode period." << endl;
                    }

                    // continue moving index forward until we reach the
                    // current macropulse time
                    do
                    {
                        index++;
                        orchardRaw[detIndex/2]->GetEntry(index);
                        procEvent.completeTime = (double)separatedEvent.extTime*pow(2,32)+separatedEvent.timetag+TIME_OFFSET;
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
            procEvent.sgQ = separatedEvent.sgQ;
            procEvent.lgQ = separatedEvent.lgQ;
            procEvent.waveform = separatedEvent.waveform;

            // only add events that come while beam is on, during the macropulse 
            double timeDiff = procEvent.completeTime-tcEvent.macroTime;

            if (timeDiff > 0 && timeDiff < MACRO_LENGTH)
            {
                orchardProcessed[detIndex/2]->Fill();
            }

            // in preparation for looping to the next event, update counters
            separatedNumberOfDPPs++;
            procEvent.evtNo++;

            prevCompleteTime = procEvent.completeTime;
            prevEvtType = separatedEvent.evtType;
            extTimePrev = separatedEvent.extTime;
            timetagPrev = separatedEvent.timetag;

            if(index%10000==0)
            {
                cout << "Event " << index << "\r";
                fflush(stdout);
            }
        }
    }
    cout << "Total number of DPP-mode events processed = " << separatedNumberOfDPPs << endl;
}

// assign target position to each waveform event
void processWaveformEvents(TFile*& sortedFile, vector<TTree*>& orchardRawW, vector<TTree*>& orchardProcessedW)
{
    TTree* targetChangerTree = (TTree*)sortedFile->Get("targetChangerTree");

    // Now that macropulses have been assigned, point macropulse variables to
    // the target changer tree
    setBranchesProcessedTC(targetChangerTree);

    for(int i=0; (size_t)i<activeWaveformChannels.size(); i++)
    {
        orchardProcessedW.push_back(new TTree((activeWaveformChannels[i]+"ProcessedTreeW").c_str(),""));
        branchProcW(orchardProcessedW[i]);
        //orchardProcessedW[i]->SetDirectory(tempFile);
    }

    long separatedNumberOfWaveforms = 0;

    for(int detIndex=0; detIndex<NUMBER_OF_CHANNELS*2; detIndex+=2)
    {
        // Attach variables to correct waveform tree
        setBranchesSeparatedW(orchardRawW[detIndex/2]);
        
        int totalEntries = orchardRawW[detIndex/2]->GetEntries();

        int prevMacroNo = 0;

        targetChangerTree->GetEntry(0);

        for(int i=0; i<totalEntries; i++)
        {
            orchardRawW[detIndex/2]->GetEntry(i);

            if(separatedEvent.evtNo==0)
            {
                // Start of new waveform mode period:
                // Move forward to the next mode change in the target changer

                do
                {
                    targetChangerTree->GetEntry(tcEvent.macroNo+1);

                } while(tcEvent.modeChange==0 && tcEvent.macroNo+1<targetChangerTree->GetEntries());

                if(tcEvent.macroNo+1>=targetChangerTree->GetEntries())
                {
                    // reached end of targetChangerTree
                    break;
                }

                prevMacroNo=tcEvent.macroNo;
            }

            // calculate this event's time
            procEvent.completeTime = (double)separatedEvent.extTime*pow(2,32)+separatedEvent.timetag;

            // assign remaining event variables to processed event
            procEvent.macroNo = tcEvent.macroNo;
            procEvent.macroTime = tcEvent.macroTime;
            procEvent.targetPos = tcEvent.targetPos;

            procEvent.evtNo = separatedEvent.evtNo;
            procEvent.waveform = separatedEvent.waveform;

            orchardProcessedW[detIndex/2]->Fill();
            separatedNumberOfWaveforms++;
        }
    }
    cout << "Total number of waveform-mode events processed = " << separatedNumberOfWaveforms << endl;
}
