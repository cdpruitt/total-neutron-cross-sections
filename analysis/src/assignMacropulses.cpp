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

#include "../include/dataStructures.h" // defines the C-structs that hold each event's data
#include "../include/branches.h" // used to map C-structs that hold raw data to ROOT trees, and vice-versa

#include "../include/assignMacropulses.h" // declarations of functions used to assign times and macropulses to events
#include "../include/experimentalConfig.h"

extern ExperimentalConfig experimentalConfig;

using namespace std;

// Use the lgQ from the target changer to determine the target position
int assignTargetPos(int lgQ)
{
    for(int i=0; (size_t)i<experimentalConfig.targetConfig.TARGET_GATES.size(); i++)
    {
        if (lgQ>=experimentalConfig.targetConfig.TARGET_GATES[i].first && lgQ<=experimentalConfig.targetConfig.TARGET_GATES[i].second)
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

void addDetectorEvent(long evtNo, TTree* detectorTree)
{
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

        SeparatedEvent separatedEvent;

        //separatedEvent.waveform = new vector<int>;
        setBranchesSeparated(rawTree, separatedEvent);

        double prevTimetag = 0;

        long totalEntries = rawTree->GetEntries();

        const double SCALEDOWN = 1;
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

        // FIX: extended time increment
        // manually increment extended time, if necessary
        /*if((prevTimetag > (double)pow(2,32) - 50000000) &&
            rawEvent.timetag < prevTimetag)
        {
           extTime++; 
        }*/

        // FIX: fine time calculation
        /*if(rawEvent.extraSelect==0)
        {
            switch(rawEvent.chNo)
            {
                case 0:
                    rawEvent.fineTime = calculateTCFineTime(
                            rawEvent.waveform,
                            rawEvent.baseline-TC_FINETIME_THRESHOLD);
                    break;
                default:
                    rawEvent.fineTime = calculateCFDTime(rawEvent.waveform,
                            rawEvent.baseline,
                            CFD_FRACTION,
                            CFD_DELAY,
                            false);
                    break;
            }
        }*/

        // FIX: assign target changer position to macropulse
        /*if(rawEvent.chNo==0)
          {
          lgQTargetChanger = rawEvent.lgQ;
          }*/

        /*if(rawEvent.chNo==1)
          {
          rawEvent.lgQ = lgQTargetChanger;
          }*/

        /*if(prevEvtType!=rawEvent.evtType)
          {
          extTime=0;
          prevTimetag=0;
          }

          rawEvent.extTime = extTime;
          */


        else if(i==1)
        {
            TFile* sortedFile = new TFile(sortedFileName.c_str(), "UPDATE");
            TTree* sortedTree = new TTree(channelMap[i].c_str(),"");

            TargetChangerEvent tcEvent;

            branchTargetChanger(sortedTree, tcEvent);

            for(int j=0; j<totalEntries; j++)
            {
                rawTree->GetEntry(j);

                tcEvent.targetPos = assignTargetPos(separatedEvent.lgQ);

                tcEvent.macroTime = experimentalConfig.timeConfig.SAMPLE_PERIOD*(pow(2,31)*separatedEvent.extTime + separatedEvent.timetag + separatedEvent.fineTime);

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

                // Check for digitizer error (incrementing extTime before clearing timetag)
                if ((separatedEvent.extTime>0) && separatedEvent.timetag > pow(2,32)-1000)
                {
                    cerr << "Found a target changer event with a timestamp-reset failure (i.e., extTime incremented before timetag was reset to 0). MacroNo = " << tcEvent.macroNo << ", extTime = " << separatedEvent.extTime << ", timetag = " << separatedEvent.timetag << endl;
                    cerr << "Skipping to next target changer event..." << endl;
                    continue;
                }

                // fill w/ new macropulse data
                tcEvent.lgQ = separatedEvent.lgQ;
                tcEvent.fineTime = separatedEvent.fineTime;

                //tcEvent.waveform = separatedEvent.waveform;

                vector<int> tempWaveform = *separatedEvent.waveform;
                tcEvent.waveform = &tempWaveform;

                sortedTree->Fill();

                // update macropulse counters
                tcEvent.macroNo++;

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
            string fineTimeHName = channelMap[i]+"fineTimeH";
            TH1I* fineTimeH = new TH1I(fineTimeHName.c_str(),fineTimeHName.c_str(),1000,-50,50);
            TFile* sortedFile = new TFile(sortedFileName.c_str(), "UPDATE");
            TTree* sortedTree = new TTree(channelMap[i].c_str(),"");
            sortedTree->SetDirectory(sortedFile);

            ProcessedEvent procEvent;
            branchProc(sortedTree, procEvent);

            TTree* targetChangerTree = (TTree*)sortedFile->Get(channelMap[1].c_str());
            TargetChangerEvent tcEvent;
            setBranchesProcessedTC(targetChangerTree, tcEvent);

            long targetChangerEntries = targetChangerTree->GetEntries();
            long currentTargetChangerEntry = 0;

            targetChangerTree->GetEntry(currentTargetChangerEntry);

            procEvent.macroNo = tcEvent.macroNo;
            procEvent.macroTime = tcEvent.macroTime;
            procEvent.targetPos = tcEvent.targetPos;
            // procEvent now has data from the first macropulse

            currentTargetChangerEntry += 1;
            targetChangerTree->GetEntry(currentTargetChangerEntry);

            long evtNo = 0;

            // The experiment setup has a cable and electronics delay (TIME_OFFSET)
            // that is unique to each channel. All times are taken relative to the
            // macropulse state time.
            double TIME_OFFSET;

            switch(i)
            {
                case 2:
                    TIME_OFFSET = experimentalConfig.timeConfig.MACROPULSE_OFFSET;
                    break;
                case 3:
                case 4:
                    TIME_OFFSET = experimentalConfig.timeConfig.MACROPULSE_OFFSET;
                    break;
                case 5:
                    TIME_OFFSET = experimentalConfig.timeConfig.MACROPULSE_OFFSET-experimentalConfig.timeConfig.VETO_OFFSET;
                    break;
                case 6:
                case 7:
                    TIME_OFFSET = experimentalConfig.timeConfig.MACROPULSE_OFFSET;
                    break;

                default:
                    cerr << "Error: non-existent channel index given by detIndex" << endl;
                    exit(1);
            }

            prevTimetag = 0;

            for(int j=0; j<totalEntries; j++)
            {
                rawTree->GetEntry(j);

                procEvent.completeTime = experimentalConfig.timeConfig.SAMPLE_PERIOD*
                    (pow(2,31)*separatedEvent.extTime +
                     separatedEvent.timetag +
                     separatedEvent.fineTime) +
                    TIME_OFFSET;

                fineTimeH->Fill(separatedEvent.fineTime);

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
                        currentTargetChangerEntry += 1;
                        targetChangerTree->GetEntry(currentTargetChangerEntry);
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
                        currentTargetChangerEntry += 1;
                        targetChangerTree->GetEntry(currentTargetChangerEntry);
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
                            currentTargetChangerEntry += 1;
                            targetChangerTree->GetEntry(currentTargetChangerEntry);
                        }

                        else
                        {
                            cout << "Reached the end of the target changer tree. Discarding all remaining events." << endl;
                            break;
                        }
                    }
                }

                // update unique event information
                procEvent.evtNo = evtNo;
                procEvent.sgQ = separatedEvent.sgQ;
                procEvent.lgQ = separatedEvent.lgQ;
                procEvent.fineTime = separatedEvent.fineTime;
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
                sortedTree->Fill();

                prevTimetag = procEvent.completeTime; 
                evtNo++;

                if(j%10000==0)
                {
                    cout << "processed " << j << " events in " << channelMap[i] << " tree.\r";
                    fflush(stdout);
                }
            }

            fineTimeH->Write();
            sortedTree->Write();
            sortedFile->Close();

        }

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

        //delete separatedEvent.waveform;
        //separatedEvent.waveform = new vector<int>;
        SeparatedEvent separatedEvent;
        setBranchesSeparated(treeToSort, separatedEvent);

        TFile* sortedFile = new TFile(sortedFileName.c_str(), "UPDATE");
        TTree* sortedTree = new TTree(treeName.c_str(),"");
        ProcessedEvent procEvent;
        branchProcW(sortedTree, procEvent);

        long totalEntries = treeToSort->GetEntries();

        TTree* targetChangerTree = (TTree*)sortedFile->Get(channelMap[1].c_str());
        TargetChangerEvent tcEvent;
        setBranchesProcessedTC(targetChangerTree, tcEvent);
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
                procEvent.completeTime = pow(2,32)*separatedEvent.extTime + experimentalConfig.timeConfig.SAMPLE_PERIOD*separatedEvent.timetag;
                procEvent.macroNo = tcEvent.macroNo;
                procEvent.macroTime = tcEvent.macroTime;
                procEvent.targetPos = tcEvent.targetPos;

                evtNo = 0;
                prevTimetag = 0;
                while(prevTimetag < procEvent.completeTime)
                {
                    if(prevTimetag+experimentalConfig.facilityConfig.MACRO_LENGTH < procEvent.completeTime)
                    {
                        evtNo=0;
                    }

                    evtNo++;
                    addDetectorEvent(evtNo, sortedTree);
                    prevTimetag = procEvent.completeTime;
                    treeToSort->GetEntry(currentWaveformEvent++);
                    procEvent.completeTime = pow(2,32)*separatedEvent.extTime + experimentalConfig.timeConfig.SAMPLE_PERIOD*separatedEvent.timetag;
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
