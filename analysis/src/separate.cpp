/******************************************************************************
                               separate.cpp 
******************************************************************************/
// separate.cpp takes a ROOT tree containing all events as input (from ./raw),
// and splits it into channel-specific ROOT trees. and assigns each event to a
// macropulse in preparation for producing cross section plots.

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

using namespace std;

extern SeparatedEvent separatedEvent;
extern TargetChangerEvent tcEvent;

long separatedNumberOfCh0Waveforms = 0;
long separatedNumberOfCh2Waveforms = 0;
long separatedNumberOfCh4Waveforms = 0;
long separatedNumberOfCh6Waveforms = 0;

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

    if(separatedEvent.evtType==2)
    {
        if(separatedEvent.chNo==0)
        {
            separatedNumberOfCh0Waveforms++;
        }

        if(separatedEvent.chNo==2)
        {
            separatedNumberOfCh2Waveforms++;
        }

        if(separatedEvent.chNo==4)
        {
            separatedNumberOfCh4Waveforms++;
        }

        if(separatedEvent.chNo==6)
        {
            separatedNumberOfCh6Waveforms++;
        }
 
    }
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
    separatedEvent.waveform->clear();

    // fill waveform with wavelet data from an event
    // (the raw tree's branch points to dummyWaveform)
    for(int k=0; (size_t)k<dummyWaveform->size(); k++)
    {
        separatedEvent.waveform->push_back(dummyWaveform->at(k));
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
    /*separatedNumberOfCh0Waveforms = 0;
    separatedNumberOfCh2Waveforms = 0;
    separatedNumberOfCh4Waveforms = 0;

    separatedNumberOfDPPs = 0;
    separatedNumberOfWaveforms = 0;
    */

    TFile* tempFile = new TFile(tempFileName.c_str(), "RECREATE");

    // Create the new empty trees
    // Each channel has a separate tree for DPP data and for waveform mode data
    for(int i=0; (size_t)i<activeDPPChannels.size(); i++)
    {
        orchardRaw.push_back(new TTree((activeDPPChannels[i]+"RawTree").c_str(),""));
        branchSplit(orchardRaw[i]);
        orchardRaw[i]->SetDirectory(tempFile);
    }

    for(int i=0; (size_t)i<activeWaveformChannels.size(); i++)
    {
        orchardRawW.push_back(new TTree((activeWaveformChannels[i]+"RawTreeW").c_str(),""));
        branchSplitW(orchardRawW[i]);
        orchardRawW[i]->SetDirectory(tempFile);
    }

    cout << "Separating events by channel and event type..." << endl;

    TFile* rawFile = new TFile(rawFileName.c_str(),"READ");
    TTree* inputTree = (TTree*)rawFile->Get("tree");

    if(!rawFile->Get("tree"))
    {
        cerr << "Error: failed to find raw tree in " << rawFileName << endl;
        exit(1);
    }

    // link the tree from the input file to our event variables
    setBranchesSeparated(inputTree);
    
    int prevEvtType = 0;
    
    int totalEntries = inputTree->GetEntries();

    totalEntries /= DEBUG_SCALEDOWN; // for debugging:
                                     // use to separate only a subset of total
                                     // events and ignore the rest

    // To uniquely identify each event, we assign each channel's events an event
    // number (evtNo), which is the event's order in its macropulse.
    vector<int> evtNo(activeDPPChannels.size());

    // Use the previous event's timetags to re-insert missing extTimes
    vector<int> currentExtTime(activeDPPChannels.size(),0);
    vector<unsigned int> prevTimetag(activeDPPChannels.size(),0);
    
    // Loop through all events in input tree and separate them into channel-
    // specific trees
    for (int index=0; index<totalEntries; index++)
    {
        inputTree->GetEntry(index);

        // if the event is a target changer event (chNo = 0),
        // processTargetChanger will handle it separately because target
        // changer events are special, non-detector events
        /*if(separatedEvent.chNo==0)
        {
            continue;
        }*/

        /*if(prevTimetag[separatedEvent.chNo/2] > 3000000000)
        {
            cout << separatedEvent.timetag << endl;
        }*/

        if(prevTimetag[separatedEvent.chNo/2] > separatedEvent.timetag
           && prevTimetag[separatedEvent.chNo/2] > pow(2,32)-50000000) // 50 ms before extTime kicks in
        {
            currentExtTime[separatedEvent.chNo/2]++;
        }

        if (separatedEvent.evtType==1)
        {
            // DPP mode

            separatedEvent.extTime = currentExtTime[separatedEvent.chNo/2];

            fillRawTree(orchardRaw[separatedEvent.chNo/2]);
        }

        if (separatedEvent.evtType==2)
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

            separatedEvent.evtNo = evtNo[separatedEvent.chNo/2];
            fillRawTree(orchardRawW[separatedEvent.chNo/2]);
            evtNo[separatedEvent.chNo/2]++;

            for(auto &value: currentExtTime)
            {
                value = 0;
            }
        }

        if(index%10000==0)
        {
            cout << "processed " << index << " events.\r";
            fflush(stdout);
        }

        prevEvtType = separatedEvent.evtType;
        prevTimetag[separatedEvent.chNo/2] = separatedEvent.timetag;
    }

    cout << "Separated " << totalEntries << " events into channels 2, 4 and 6." << endl;
    tempFile->Write();

    rawFile->Close();
}

// assign a macropulse, target position, etc to every target changer event
void processTargetChanger(string rawFileName, TFile*& sortedFile)
{
    TTree* targetChangerTree = new TTree("targetChangerTree","");
    branchTargetChanger(targetChangerTree);
    targetChangerTree->SetDirectory(sortedFile);

    TFile* rawData = new TFile(rawFileName.c_str(),"READ");
    TTree* inputTree = (TTree*)rawData->Get("tree");
    setBranchesSeparated(inputTree);

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

    int currentExtTime = 0;

    for(int i=0; i<totalEntries; i++)
    {
        inputTree->GetEntry(i);

        if(separatedEvent.chNo==0)
        {
            // treat DPP and waveform mode events differently
            if(separatedEvent.evtType==1)
            {
                // DPP mode
                if (separatedEvent.lgQ==65535)
                {
                    // ignore target changer events having an off-scale integrated
                    // charge - these are suspected to be retriggers and are NOT
                    // time-correlated with the facility RF time reference.
                    continue;
                }

                if(timetagPrev > separatedEvent.timetag
                && timetagPrev > pow(2,32)-50000000) // 50 ms before extTime kicks in
                {
                    currentExtTime++;
                }

                separatedEvent.extTime = currentExtTime;

                // assign the macropulse start time
                tcEvent.macroTime = (double)separatedEvent.extTime*pow(2,32)+separatedEvent.timetag;

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
                if (separatedEvent.extTime > extTimePrev && separatedEvent.timetag > pow(2,32)-1000)
                {
                    cerr << "Found a target changer event with a timestamp-reset failure (i.e., extTime incremented before timetag was reset to 0). MacroNo = " << tcEvent.macroNo << ", extTime = " << separatedEvent.extTime << ", extTimePrev = " << extTimePrev << ", timetag = " << separatedEvent.timetag << ", timetagPrev = " << timetagPrev << endl;
                    cerr << "Skipping to next target changer event..." << endl;
                    continue;
                }

                tcEvent.targetPos = assignTargetPos(separatedEvent.lgQ);

                /*************************************************************/
                // all the variables are updated: fill tree
                targetChangerTree->Fill();
                tcEvent.macroNo++;

                // before we pull the next event, save the previous event's data so we can
                // compare them with the new event's data
                extTimePrev = separatedEvent.extTime;
                timetagPrev = separatedEvent.timetag;

                if(tcEvent.macroNo%100==0)
                {
                    cout << "Target position = " << tcEvent.targetPos << ", macroNo = " << tcEvent.macroNo << "\r";
                    fflush(stdout);
                }
            }

            else
            {
                currentExtTime = 0;
            }

            prevEvtType = separatedEvent.evtType;
            // move to next event in the loop

            //cerr << "macroNo = " << tcEvent.macroNo << ", macroTime = " << tcEvent.macroTime
            //    << ", modeChange = " << tcEvent.modeChange << endl;
        }
    }

    // loop finished
    // all target changer events separated into targetChanger and targetChangerTreeW
    // time to move to detector channels

    cout << endl << "... done." << endl;

    sortedFile->cd();
    targetChangerTree->Write();
}
