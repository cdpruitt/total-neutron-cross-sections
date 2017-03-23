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
#include <utility>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#include "../include/dataStructures.h"
#include "../include/analysisConstants.h"
#include "../include/physicalConstants.h"
#include "../include/separate.h"
#include "../include/branches.h"

using namespace std;

extern SeparatedEvent separatedEvent;
extern ProcessedEvent procEvent;
extern TargetChangerEvent tcEvent;

void fillRawTreeW(TTree* tree)
{
    // To save only 1 DPP event's waveform in every 10000 events, uncomment this
    // (Note - will reduce output tree size by ~half)
    /*if ((int)completeTime%10000 != 0 && evtType==1)
      {
      waveform->clear(); 
      }*/

    tree->Fill();
}

// Use the lgQ from the target changer to determine the target position
int assignTargetPos(int lgQ)
{
    for(int i=0; (size_t)i<tarGates.size(); i++)
    {
        if (lgQ>tarGates[i].first && lgQ<tarGates[i].second)
        {
            // lgQ fits within this gate
            return i; // target positions start from 1
        }
    }

    // lgQ doesn't fit within one of the lgQ gates, so set
    // the target position to 0 (discarded in later analysis).
    cerr << "Error: lgQ outside target positions." << endl;
    exit(1);
}

void addTCEvent(TTree* targetChangerTree)
{
    // fill w/ new macropulse data
    tcEvent.targetPos = assignTargetPos(separatedEvent.lgQ);
    tcEvent.lgQ = separatedEvent.lgQ;

    targetChangerTree->Fill();

    // update macropulse counters
    tcEvent.macroNo++;
}

void addDetectorEvent(long evtNo, TTree* detectorTree)
{
    // update unique event information
    procEvent.evtNo = evtNo;
    procEvent.sgQ = separatedEvent.sgQ;
    procEvent.lgQ = separatedEvent.lgQ;
    procEvent.waveform = separatedEvent.waveform;

    // only add events that come while beam is on, during the macropulse 
/*    double timeDiff = procEvent.completeTime-tcEvent.macroTime;
    if (timeDiff < 0 && timeDiff > MACRO_LENGTH)
    {
        return;
    }
    */

    // fill tree w/ event data

    /*if(tcEvent.macroNo > 12090 && tcEvent.macroNo < 12095)
    {
        cout << "macroNo = " << tcEvent.macroNo << ", macroTime = " << tcEvent.macroTime
            << ", evtNo = " << procEvent.evtNo << ", channel = " << chNo <<  ", completeTime = " << procEvent.completeTime << endl;
    }*/
    detectorTree->Fill();
}

void assignMacropulses(string rawFileName, string sortedFileName, vector<string> channelMap)
{
    TFile* rawFile = new TFile(rawFileName.c_str(),"READ");
    if (!rawFile)
    {
        cerr << "Failed to open " << rawFileName << ". Please check that the file exists" << endl;
        exit(1);
    }
    cout << rawFileName << " opened successfully. Start reading events..." << endl;

    TFile* sortedFile = new TFile(sortedFileName.c_str(), "RECREATE");

    TH1I* fineTimeHBlank = new TH1I("fineTimeHBlank","fineTimeHBlank",6000,-3,3);
    TH1I* fineTimeHTarget1 = new TH1I("fineTimeHTarget1","fineTimeHTarget1",6000,-3,3);
    TH1I* fineTimeHTarget2 = new TH1I("fineTimeHTarget2","fineTimeHTarget2",6000,-3,3);
    TH1I* fineTimeHTarget3 = new TH1I("fineTimeHTarget3","fineTimeHTarget3",6000,-3,3);
    TH1I* fineTimeHTarget4 = new TH1I("fineTimeHTarget4","fineTimeHTarget4",6000,-3,3);
    TH1I* fineTimeHTarget5 = new TH1I("fineTimeHTarget5","fineTimeHTarget5",6000,-3,3);

    for(int i=0; i<channelMap.size(); i++)
    {
        if(channelMap[i]=="-")
        {
            continue;
        }

        TTree* rawTree = (TTree*)rawFile->Get(channelMap[i].c_str());
        if(!rawTree)
        {
            cerr << "Error: couldn't find channel " << i << " tree when attempting to assign macropulses. Exiting... " << endl;
            exit(1);
        }

        setBranchesSeparated(rawTree);

        TTree* sortedTree = new TTree(channelMap[i].c_str(),"");
        sortedTree->SetDirectory(sortedFile);

        if(i==0)
        {
            // first channel is the target changer
            branchTargetChanger(sortedTree);
        }

        else
        {
            branchProc(sortedTree);
        }

        double prevTimetag = 0;

        long totalEntries = rawTree->GetEntries();

        totalEntries /= SCALEDOWN; // for debugging:
        // use to separate only a subset of total
        // events and ignore the rest

        if(i==0)
        {
            for(int j=0; j<totalEntries; j++)
            {
                rawTree->GetEntry(j);

                tcEvent.macroTime = (double)pow(2,32)*separatedEvent.extTime + separatedEvent.timetag;

                if(prevTimetag > tcEvent.macroTime)
                {
                    // mode change
                    tcEvent.modeChange = 1;
                }

                else
                {
                    tcEvent.modeChange = 0;
                }

                addTCEvent(sortedTree);

                prevTimetag = tcEvent.macroTime;

                if(j%1000==0)
                {
                    cout << "processed " << j << " events in target changer tree.\r";
                    fflush(stdout);
                }

            }
        }

        else
        {
            TTree* targetChangerTree = (TTree*)sortedFile->Get(channelMap[0].c_str());
            setBranchesProcessedTC(targetChangerTree);
            long targetChangerEntries = targetChangerTree->GetEntries();
            long currentTargetChangerEntry = 0;
            targetChangerTree->GetEntry(currentTargetChangerEntry);

            procEvent.macroNo = tcEvent.macroNo;
            procEvent.macroTime = tcEvent.macroTime;
            procEvent.targetPos = tcEvent.targetPos;
            // procEvent now has data from the first macropulse

            targetChangerTree->GetEntry(++currentTargetChangerEntry);

            long evtNo = 0;

            // The experiment setup has a cable and electronics delay (TIME_OFFSET)
            // that is unique to each channel. All times are taken relative to the
            // macropulse state time.
            double TIME_OFFSET;

            switch(i)
            {
                case 2:
                    TIME_OFFSET = MACROPULSE_OFFSET;
                    break;
                case 3:
                case 4:
                case 5:
                    TIME_OFFSET = MACROPULSE_OFFSET;
                    break;
                case 6:
                    TIME_OFFSET = MACROPULSE_OFFSET-VETO_OFFSET;
                    break;
                default:
                    cerr << "Error: non-existent channel index given by detIndex" << endl;
                    exit(1);
            }

            for(int j=0; j<totalEntries; j++)
            {
                rawTree->GetEntry(j);

                procEvent.completeTime = (double)pow(2,32)*separatedEvent.extTime
                    + separatedEvent.timetag + TIME_OFFSET;

                if(i==4 || i==5)
                {
                    procEvent.completeTime += separatedEvent.fineTime;

                    switch(procEvent.targetPos)
                    {
                        case 1:
                            fineTimeHBlank->Fill(separatedEvent.fineTime);
                            break;
                        case 2:
                            fineTimeHTarget1->Fill(separatedEvent.fineTime);
                            break;
                        case 3:
                            fineTimeHTarget2->Fill(separatedEvent.fineTime);
                            break;
                        case 4:
                            fineTimeHTarget3->Fill(separatedEvent.fineTime);
                            break;
                        case 5:
                            fineTimeHTarget4->Fill(separatedEvent.fineTime);
                            break;
                        case 6:
                            fineTimeHTarget5->Fill(separatedEvent.fineTime);
                            break;
                        default:
                            break;
                    }
                }

                    if(procEvent.completeTime < prevTimetag)
                {
                    // Cycle the macropulse forward to the next DPP/waveform
                    // mode change, if necessary

                    while(currentTargetChangerEntry+1<targetChangerEntries
                            && tcEvent.macroTime > procEvent.macroTime) 
                    {
                        procEvent.macroNo = tcEvent.macroNo;
                        procEvent.macroTime = tcEvent.macroTime;
                        procEvent.targetPos = tcEvent.targetPos;
                        targetChangerTree->GetEntry(++currentTargetChangerEntry);
                    }

                    /*                    // DPP/waveform mode change
                                          procEvent.macroNo = tcEvent.macroNo;
                                          procEvent.macroTime = tcEvent.macroTime;
                                          procEvent.targetPos = tcEvent.targetPos;
                                          evtNo = 0;
                                          prevTimetag = 0;

                                          if(currentTargetChangerEntry+1<targetChangerEntries)
                                          {

                                          targetChangerTree->GetEntry(++currentTargetChangerEntry);
                                          }

                                          else
                                          {
                                          cout << "Reached the end of the target changer tree. Discarding all remaining events." << endl;
                                          break;
                                          }
                                          */

                    // macropulse is now the first in next DPP/waveform change

                    while(currentTargetChangerEntry+1<targetChangerEntries &&
                            tcEvent.macroTime < procEvent.macroTime)
                    {
                        procEvent.macroNo = tcEvent.macroNo;
                        procEvent.macroTime = tcEvent.macroTime;
                        procEvent.targetPos = tcEvent.targetPos;
                        targetChangerTree->GetEntry(++currentTargetChangerEntry);
                    }
                }

                if(procEvent.completeTime > tcEvent.macroTime)
                {
                    if(tcEvent.macroTime > procEvent.macroTime)
                    {
                        // macropulse change
                        procEvent.macroNo = tcEvent.macroNo;
                        procEvent.macroTime = tcEvent.macroTime;
                        procEvent.targetPos = tcEvent.targetPos;
                        evtNo = 0;

                        if(currentTargetChangerEntry+1<targetChangerEntries)
                        {
                            targetChangerTree->GetEntry(++currentTargetChangerEntry);
                        }

                        else
                        {
                            cout << "Reached the end of the target changer tree. Discarding all remaining events." << endl;
                            break;
                        }
                    }
                }

                addDetectorEvent(evtNo, sortedTree);

                prevTimetag = procEvent.completeTime; 
                evtNo++;

                if(j%10000==0)
                {
                    cout << "processed " << j << " events in " << channelMap[i] << " tree.\r";
                    fflush(stdout);
                }
            }
        }

        // Check for digitizer error (incrementing extTime before clearing timetag)
        /*if (extTime[separatedEvent.chNo] > extTimePrev && separatedEvent.timetag > pow(2,32)-1000)
          {
          cerr << "Found a target changer event with a timestamp-reset failure (i.e., extTime incremented before timetag was reset to 0). MacroNo = " << tcEvent.macroNo << ", extTime = " << separatedEvent.extTime << ", extTimePrev = " << extTimePrev << ", timetag = " << separatedEvent.timetag << ", timetagPrev = " << timetagPrev << endl;
          cerr << "Skipping to next target changer event..." << endl;
          continue;
          }*/

        sortedTree->Write();

        cout << endl;
    }

    for(int i=0; i<channelMap.size(); i++)
    {
        if(channelMap[i]=="-")
        {
            continue;
        }

        string treeName = channelMap[i]+"W";
        TTree* treeToSort = (TTree*)rawFile->Get(treeName.c_str());
        if(!treeToSort)
        {
            cerr << "Error: couldn't find channel " << i << " tree when attempting to assign macropulses. Exiting... " << endl;
            exit(1);
        }

        setBranchesSeparated(treeToSort);

        TTree* sortedTree = new TTree(treeName.c_str(),"");
        sortedTree->SetDirectory(sortedFile);
        branchProcW(sortedTree);

        long totalEntries = treeToSort->GetEntries();

        TTree* targetChangerTree = (TTree*)sortedFile->Get(channelMap[0].c_str());
        setBranchesProcessedTC(targetChangerTree);
        long targetChangerEntries = targetChangerTree->GetEntries();

        // procEvent now has data from the first macropulse

        long evtNo = 0;
        long currentWaveformEvent = -1;
        double prevTimetag = 0;

        for(int j=0; j<targetChangerEntries; j++)
        {
            targetChangerTree->GetEntry(j);

            if(tcEvent.modeChange==1)
            {
                treeToSort->GetEntry(currentWaveformEvent++);
                procEvent.completeTime = pow(2,32)*separatedEvent.extTime + separatedEvent.timetag;
                procEvent.macroNo = tcEvent.macroNo;
                procEvent.macroTime = tcEvent.macroTime;
                procEvent.targetPos = tcEvent.targetPos;

                evtNo = 0;
                prevTimetag = 0;
                while(prevTimetag < procEvent.completeTime)
                {
                    if(prevTimetag+MACRO_LENGTH < procEvent.completeTime)
                    {
                        evtNo=0;
                    }

                    evtNo++;
                    addDetectorEvent(evtNo, sortedTree);
                    prevTimetag = procEvent.completeTime;
                    treeToSort->GetEntry(currentWaveformEvent++);
                    procEvent.completeTime = pow(2,32)*separatedEvent.extTime + separatedEvent.timetag;
                }
            }

            if(currentWaveformEvent%10==0)
            {
                cout << "processed " << j << "waveform events in " << channelMap[i] << " tree.\r";
                fflush(stdout);

            }
        }

        sortedTree->Write();
    }

    cout << endl;

    fineTimeHBlank->Write();
    fineTimeHTarget1->Write();
    fineTimeHTarget2->Write();
    fineTimeHTarget3->Write();
    fineTimeHTarget4->Write();
    fineTimeHTarget5->Write();

    sortedFile->Close();
}

// Populate events from the input tree into channel-specific trees.
void separateByChannel(string rawFileName, string sortedFileName, vector<string>& channelMap)
{
    cout << "Separating events by channel and event type..." << endl;

    // open the input file
    TFile* rawFile = new TFile(rawFileName.c_str(),"READ");
    TTree* inputTree = (TTree*)rawFile->Get("tree");

    if(!rawFile->Get("tree"))
    {
        cerr << "Error: failed to find raw tree in " << rawFileName << endl;
        exit(1);
    }

    // if no channel mapping can be established for this run, no separation can
    // be done - exit
    if(!channelMap.size())
    {
        cerr << "Error: failed to establish channel mapping for this run. Exiting..." << endl;
        exit(1);
    }

    // link the tree from the input file to our event variables
    setBranchesSeparated(inputTree);

    int totalEntries = inputTree->GetEntries();

    totalEntries /= SCALEDOWN; // for debugging:
                               // use to separate only a subset of total
                               // events and ignore the rest

    // create an output file
    TFile* sortedFile = new TFile(sortedFileName.c_str(),"CREATE");
    vector<TTree*> orchardProcessed; // channel-specific DPP events assigned to macropulses
    vector<TTree*> orchardProcessedW;// channel-specific waveform events assigned to macropulses

    // Create trees to be filled with sorted data
    // Each channel has a separate tree for DPP data and for waveform mode data
    for(int i=0; (size_t)i<channelMap.size(); i++)
    {
        orchardProcessed.push_back(new TTree((channelMap[i]).c_str(),""));
        orchardProcessedW.push_back(new TTree((channelMap[i]+"W").c_str(),""));

        if(i==0)
        {
            branchTargetChanger(orchardProcessed[i]);
            orchardProcessed[i]->SetDirectory(sortedFile);
            continue;
        }

        if(channelMap[i]=="-")
        {
            continue;
        }

        branchProc(orchardProcessed[i]);
        orchardProcessed[i]->SetDirectory(sortedFile);

        branchProcW(orchardProcessedW[i]);
        orchardProcessedW[i]->SetDirectory(sortedFile);
    }

    for (int index=0; index<totalEntries; index++)
    {
        inputTree->GetEntry(index);
        
        if(separatedEvent.evtType==1)
        {
            // DPP mode
            orchardProcessed[separatedEvent.chNo]->Fill();

            /*case 0:
            // target changer event
            addTCEvent(evtNo, extTime, orchardProcessed[separatedEvent.chNo]);
            break;
            case 2:
            // clear finetime for monitor events!
            separatedEvent.fineTime=0;
            case 3:
            case 4:
            case 5:
            case 6:
            addDetectorEvent(evtNo, extTime, separatedEvent.chNo, 
            orchardProcessed[separatedEvent.chNo]);
            break;
            */
        }

        else if(separatedEvent.evtType==2)
        {
            cout << separatedEvent.evtType << endl;
            // Waveform mode
            orchardProcessedW[separatedEvent.chNo]->Fill();
        }
    }

    cout << "Separated " << totalEntries << " events into channels." << endl;

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

    sortedFile->Close();
}
