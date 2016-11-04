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
extern ProcessedEvent procEvent;
extern TargetChangerEvent tcEvent;

long separatedNumberOfCh0Waveforms = 0;
long separatedNumberOfCh2Waveforms = 0;
long separatedNumberOfCh4Waveforms = 0;
long separatedNumberOfCh6Waveforms = 0;

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

void addTCEvent(vector<int>& evtNo, vector<int>& extTime, TTree* targetChangerTree)
{
    // reset event counter for each channel
    for(auto &value: evtNo)
    {
        value = 0;
    }

    // fill w/ new macropulse data
    tcEvent.macroTime = pow(2,32)*extTime[0] + separatedEvent.timetag;
    tcEvent.targetPos = assignTargetPos(separatedEvent.lgQ);

    targetChangerTree->Fill();

    // update macropulse counters
    tcEvent.macroNo++;
    tcEvent.modeChange = 0;

    if(tcEvent.macroNo%100==0)
    {
        cout << "Target position = " << tcEvent.targetPos << ", macroNo = " << tcEvent.macroNo << "\r";
        fflush(stdout);
    }
}

void addDetectorEvent(vector<int>& evtNo, vector<int>& extTime, int chNo, TTree* detectorTree)
{
    // The experiment setup has a cable and electronics delay (TIME_OFFSET)
    // that is unique to each channel. All times are taken relative to the
    // macropulse state time.
    double TIME_OFFSET;

    switch(chNo)
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


    // update macropulse information
    procEvent.macroNo = tcEvent.macroNo;
    procEvent.macroTime = tcEvent.macroTime;
    procEvent.targetPos = tcEvent.targetPos;

    // update unique event information
    procEvent.evtNo = evtNo[chNo/2];
    procEvent.completeTime = pow(2,32)*extTime[chNo/2] + separatedEvent.timetag + separatedEvent.fineTime + TIME_OFFSET;
    procEvent.sgQ = separatedEvent.sgQ;
    procEvent.lgQ = separatedEvent.lgQ;
    procEvent.waveform = separatedEvent.waveform;

    // only add events that come while beam is on, during the macropulse 
    double timeDiff = procEvent.completeTime-tcEvent.macroTime;
    if (timeDiff < 0 && timeDiff > MACRO_LENGTH)
    {
        return;
    }

    // fill tree w/ event data
    detectorTree->Fill();
}

// Populate events from the input tree into channel-specific trees.
void separateByChannel(string rawFileName, string sortedFileName, vector<TTree*>& orchardProcessed, vector<TTree*>& orchardProcessedW)
{
    /*separatedNumberOfCh0Waveforms = 0;
    separatedNumberOfCh2Waveforms = 0;
    separatedNumberOfCh4Waveforms = 0;

    separatedNumberOfDPPs = 0;
    separatedNumberOfWaveforms = 0;
    */

    TFile* sortedFile = new TFile(sortedFileName.c_str(), "RECREATE");

    // Create the new empty trees
    // Each channel has a separate tree for DPP data and for waveform mode data
    for(int i=0; (size_t)i<activeDPPChannels.size(); i++)
    {
        orchardProcessed.push_back(new TTree((activeDPPChannels[i]+"ProcessedTree").c_str(),""));
        if(i==0)
        {
            branchTargetChanger(orchardProcessed[i]);
            orchardProcessed[i]->SetDirectory(sortedFile);
            continue;
        }
        branchProc(orchardProcessed[i]);
        orchardProcessed[i]->SetDirectory(sortedFile);
    }

    for(int i=0; (size_t)i<activeWaveformChannels.size(); i++)
    {
        orchardProcessedW.push_back(new TTree((activeWaveformChannels[i]+"ProcessedTreeW").c_str(),""));
        branchProcW(orchardProcessedW[i]);
        orchardProcessedW[i]->SetDirectory(sortedFile);
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
    
    int totalEntries = inputTree->GetEntries();

    totalEntries /= SCALEDOWN; // for debugging:
                                     // use to separate only a subset of total
                                     // events and ignore the rest

    // To uniquely identify each event, we assign each channel's events an event
    // number (evtNo), which is the event's order in its macropulse.
    vector<int> evtNo(activeDPPChannels.size(),0);

    // keep track of event type for each channel
    vector<int> prevEvtType(activeDPPChannels.size(),1);

    // Use the previous event's timetags to re-insert missing extTimes
    vector<int> extTime(activeDPPChannels.size(),0);
    vector<unsigned int> prevTimetag(activeDPPChannels.size(),0);
    
    // Loop through all events in input tree and separate them into channel-
    // specific trees
    for (int index=0; index<totalEntries; index++)
    {
        inputTree->GetEntry(index);

        if(prevTimetag[separatedEvent.chNo/2] > separatedEvent.timetag &&
           prevTimetag[separatedEvent.chNo/2] > pow(2,32)-20000000) // 20 ms before extTime kicks in
        {
            extTime[separatedEvent.chNo/2]++;
        }

        // Check for digitizer error (incrementing extTime before clearing timetag)
        /*if (extTime[separatedEvent.chNo/2] > extTimePrev && separatedEvent.timetag > pow(2,32)-1000)
        {
            cerr << "Found a target changer event with a timestamp-reset failure (i.e., extTime incremented before timetag was reset to 0). MacroNo = " << tcEvent.macroNo << ", extTime = " << separatedEvent.extTime << ", extTimePrev = " << extTimePrev << ", timetag = " << separatedEvent.timetag << ", timetagPrev = " << timetagPrev << endl;
            cerr << "Skipping to next target changer event..." << endl;
            continue;
        }*/

        if(prevEvtType[separatedEvent.chNo/2]!=separatedEvent.evtType)
        {
            // mode change
            evtNo[separatedEvent.chNo/2] = 0;
            extTime[separatedEvent.chNo/2] = 0;
            prevTimetag[separatedEvent.chNo/2] = 0;
            tcEvent.modeChange = 1;
        }

        if(separatedEvent.evtType==1)
        {
            // DPP mode
            switch(separatedEvent.chNo)
            {
                case 0:
                    // target changer event
                    addTCEvent(evtNo, extTime, orchardProcessed[0]);
                    break;
                case 2:
                case 4:
                case 6:
                    addDetectorEvent(evtNo, extTime, separatedEvent.chNo, 
                            orchardProcessed[separatedEvent.chNo/2]);
                    // main detector event
                    break;
                default:
                    cout << "Error: could not process event on channel " << separatedEvent.chNo << endl;
            }
        }

        else if(separatedEvent.evtType==2)
        {
            // Waveform mode
            fillRawTree(orchardProcessedW[separatedEvent.chNo/2]);
        }

        prevEvtType[separatedEvent.chNo/2] = separatedEvent.evtType;
        prevTimetag[separatedEvent.chNo/2] = procEvent.completeTime; 
        evtNo[separatedEvent.chNo/2]++;

        if(index%10000==0)
        {
            cout << "processed " << index << " events.\r";
            fflush(stdout);
        }
    }

    cout << "Separated " << totalEntries << " events into channels 2, 4 and 6." << endl;

    rawFile->Close();
    sortedFile->cd();

    for(TTree* tree : orchardProcessed)
    {
        tree->Write();
    }

    for(TTree* tree : orchardProcessedW)
    {
        tree->Write();
    }

}
