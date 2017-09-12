/******************************************************************************
  assignMacropulses.cpp 
 ******************************************************************************/
// This method takes a ROOT tree containing all events as input (from ./raw),
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
#include "../include/assignMacropulses.h"
#include "../include/branches.h"

using namespace std;

extern SeparatedEvent separatedEvent;
extern ProcessedEvent procEvent;
extern TargetChangerEvent tcEvent;

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
    cerr << endl;
    cerr << "Error: lgQ " << lgQ << " outside target positions." << endl;
    return 0;
}

void addTCEvent(TTree* targetChangerTree)
{
    // fill w/ new macropulse data
    tcEvent.lgQ = separatedEvent.lgQ;
    tcEvent.fineTime = separatedEvent.fineTime;

    //tcEvent.waveform = separatedEvent.waveform;

    vector<int> tempWaveform = *separatedEvent.waveform;
    tcEvent.waveform = &tempWaveform;

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
    /*double timeDiff = procEvent.completeTime-tcEvent.macroTime;
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

        //delete tcEvent.waveform;
        //tcEvent.waveform = new vector<int>;

        //separatedEvent.waveform = new vector<int>;
        setBranchesSeparated(rawTree);

        double prevTimetag = 0;

        long totalEntries = rawTree->GetEntries();

        totalEntries /= SCALEDOWN; // for debugging:
        // use to separate only a subset of total
        // events and ignore the rest

        if(i==0)
        {
            TFile* sortedFile = new TFile(sortedFileName.c_str(), "RECREATE");
            TTree* sortedTree = new TTree(channelMap[i].c_str(),"");
            sortedTree->Write();
            sortedFile->Close();
            continue;
        }

        /*if(separatedEvent.waveform->size() > 48)
          {
          cerr << "separatedEvent.waveform size > 48. i = " << i << endl;
          exit(1);
          }*/

        else if(i==1)
        {
            TFile* sortedFile = new TFile(sortedFileName.c_str(), "UPDATE");
            TTree* sortedTree = new TTree(channelMap[i].c_str(),"");

            branchTargetChanger(sortedTree);

            for(int j=0; j<totalEntries; j++)
            {
                rawTree->GetEntry(j);

                tcEvent.targetPos = assignTargetPos(separatedEvent.lgQ);

                tcEvent.macroTime = SAMPLE_PERIOD*(pow(2,31)*separatedEvent.extTime + separatedEvent.timetag + separatedEvent.fineTime);

                /*if(tcEvent.targetPos > 0)
                  {
                  tcEvent.macroTime += MACROTIME_TARGET_DRIFT[tcEvent.targetPos-1];
                  }*/

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

            sortedTree->Write();
            sortedFile->Close();
        }

        else
        {
            TFile* sortedFile = new TFile(sortedFileName.c_str(), "UPDATE");
            TTree* sortedTree = new TTree(channelMap[i].c_str(),"");
            sortedTree->SetDirectory(sortedFile);

            branchProc(sortedTree);

            TTree* targetChangerTree = (TTree*)sortedFile->Get(channelMap[1].c_str());
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
                    TIME_OFFSET = MACROPULSE_OFFSET;
                    break;
                case 5:
                    TIME_OFFSET = MACROPULSE_OFFSET-VETO_OFFSET;
                    break;
                case 6:
                case 7:
                    TIME_OFFSET = MACROPULSE_OFFSET;
                    break;

                default:
                    cerr << "Error: non-existent channel index given by detIndex" << endl;
                    exit(1);
            }

            prevTimetag = 0;

            for(int j=0; j<totalEntries; j++)
            {
                rawTree->GetEntry(j);

                procEvent.completeTime = SAMPLE_PERIOD*
                    (pow(2,31)*separatedEvent.extTime +
                     separatedEvent.timetag +
                     separatedEvent.fineTime) +
                    TIME_OFFSET;

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

            sortedTree->Write();
            sortedFile->Close();

        }

        // Check for digitizer error (incrementing extTime before clearing timetag)
        /*if (extTime[separatedEvent.chNo] > extTimePrev && separatedEvent.timetag > pow(2,32)-1000)
          {
          cerr << "Found a target changer event with a timestamp-reset failure (i.e., extTime incremented before timetag was reset to 0). MacroNo = " << tcEvent.macroNo << ", extTime = " << separatedEvent.extTime << ", extTimePrev = " << extTimePrev << ", timetag = " << separatedEvent.timetag << ", timetagPrev = " << timetagPrev << endl;
          cerr << "Skipping to next target changer event..." << endl;
          continue;
          }*/

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

        delete separatedEvent.waveform;
        separatedEvent.waveform = new vector<int>;
        setBranchesSeparated(treeToSort);

        TFile* sortedFile = new TFile(sortedFileName.c_str(), "UPDATE");
        TTree* sortedTree = new TTree(treeName.c_str(),"");
        branchProcW(sortedTree);

        long totalEntries = treeToSort->GetEntries();

        TTree* targetChangerTree = (TTree*)sortedFile->Get(channelMap[1].c_str());
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
                procEvent.completeTime = pow(2,32)*separatedEvent.extTime + SAMPLE_PERIOD*separatedEvent.timetag;
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
                    procEvent.completeTime = pow(2,32)*separatedEvent.extTime + SAMPLE_PERIOD*separatedEvent.timetag;
                }
            }

            if(currentWaveformEvent%10==0)
            {
                cout << "processed " << j << "waveform events in " << channelMap[i] << " tree.\r";
                fflush(stdout);

            }
        }

        sortedTree->Write();
        sortedFile->Close();
    }
}
