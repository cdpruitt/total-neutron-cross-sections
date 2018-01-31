/******************************************************************************
  assignEventsToMacropulses.cpp 
 ******************************************************************************/
// This method takes a ROOT tree containing all events as input (from ./raw),
// . and assigns each event to a
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

#include "../include/assignEventsToMacropulses.h" // declarations of functions used to assign times and macropulses to events
#include "../include/config.h"

using namespace std;

extern Config config;

int assignEventsToMacropulses(string inputFileName, string outputFileName, ofstream& logFile, vector<MacropulseEvent>& macropulseList)
{
    /**************************************************************************/
    // open input detector tree
    TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
    if(!inputFile)
    {
        cerr << "Failed to open " << inputFileName << ". Please check that the file exists" << endl;
        exit(1);
    }

    TTree* inputTree = (TTree*)inputFile->Get(config.analysis.DPP_TREE_NAME.c_str());
    if(!inputTree)
    {
        cerr << "Error: couldn't find detector tree " << config.analysis.DPP_TREE_NAME << " in "
             << inputFileName << " when attempting to assign events to macropulses. Exiting... " << endl;
        return 1;
    }
    /**************************************************************************/

    if(macropulseList.size()==0)
    {
        cerr << "Error: macropulse list size was zero in assignEventsToMacropulses. Exiting..." << endl;
        return 1;
    }

    /**************************************************************************/
    // create struct for holding event data and link it to trees
    DetectorEvent detectorEvent;
    
    unsigned int chNo;
    unsigned int cycleNumber;
    unsigned int sgQ;
    unsigned int lgQ;

    vector<int>* waveformPointer = 0;

    inputTree->SetBranchAddress("chNo",&chNo);
    inputTree->SetBranchAddress("cycleNumber",&cycleNumber);
    inputTree->SetBranchAddress("completeTime",&detectorEvent.completeTime);
    inputTree->SetBranchAddress("fineTime",&detectorEvent.fineTime);
    inputTree->SetBranchAddress("sgQ",&sgQ);
    inputTree->SetBranchAddress("lgQ",&lgQ);
    inputTree->SetBranchAddress("waveform",&waveformPointer);

    /**************************************************************************/

    vector<vector<DetectorEvent>> allEvents;

    for(auto& channel : config.digitizer.CHANNEL_MAP)
    {
        allEvents.push_back(vector<DetectorEvent>());
    }

    int inputTreeEntries = inputTree->GetEntries();
    int currentTreeEntry = 0;

    while(currentTreeEntry<inputTreeEntries)
    {
        inputTree->GetEntry(currentTreeEntry);

        DetectorEvent de = DetectorEvent(detectorEvent);

        // assign unsigned int variables -> int variables
        de.cycleNumber = cycleNumber;
        de.sgQ = sgQ;
        de.lgQ = lgQ;

        de.waveform = *waveformPointer;
        allEvents[chNo].push_back(de);
        detectorEvent.eventNo++;

        if(currentTreeEntry%10000==0)
        {
            cout << "Read " << currentTreeEntry << " \"" << config.analysis.DPP_TREE_NAME << "\" events...\r";
            fflush(stdout);
        }

        currentTreeEntry++;
    }

    TFile* outputFile = new TFile(outputFileName.c_str(),"UPDATE");

    for(auto& channel : config.digitizer.CHANNEL_MAP)
    {
        if(
                channel.second == "-" ||
                channel.second == "macroTime" ||
                channel.second == "targetChanger"
          )
        {
            continue;
        }

        cout << endl << "Start assigning \"" << channel.second
            << "\" events to macropulses..." << endl;

        /*********************************************************************/
        // create a new tree for holding sorted detector events
        TTree* outputTree = new TTree(channel.second.c_str(),"");
        if(!outputTree)
        {
            cerr << "Error: couldn't create detector tree " << channel.second
                << " when attempting to assign events to macropulses. Exiting... " << endl;
            return 1;
        }
        /*********************************************************************/

        // prepare output tree for filling
        outputTree->Branch("macroTime",&detectorEvent.macroTime,"macroTime/D");
        outputTree->Branch("completeTime",&detectorEvent.completeTime,"completeTime/d");
        outputTree->Branch("fineTime",&detectorEvent.fineTime,"fineTime/d");

        outputTree->Branch("macroNo",&detectorEvent.macroNo,"macroNo/I");
        outputTree->Branch("cycleNumber",&detectorEvent.cycleNumber,"cycleNumber/I");
        outputTree->Branch("targetPos",&detectorEvent.targetPos,"targetPos/I");
        outputTree->Branch("eventNo",&detectorEvent.eventNo,"eventNo/I");
        outputTree->Branch("sgQ",&detectorEvent.sgQ,"sgQ/I");
        outputTree->Branch("lgQ",&detectorEvent.lgQ,"lgQ/I");
        outputTree->Branch("waveform",&detectorEvent.waveform);

        int currentEvent = 0;
        int currentMacropulse = 0;

        bool endAssignment = false;

        int maxCycleNumber = allEvents[channel.first].back().cycleNumber;

        while(allEvents[channel.first][currentEvent].cycleNumber<=maxCycleNumber)
        {
            // if the macropulse's cycle number is behind, move to the next macropulse
            while(macropulseList[currentMacropulse].cycleNumber <
                    allEvents[channel.first][currentEvent].cycleNumber)
            {
                currentMacropulse++;
                if(currentMacropulse>=macropulseList.size())
                {
                    cout << "Reached end of macrotime list; end event assignment." << endl;
                    endAssignment = true;
                    break;
                }

                detectorEvent.eventNo=0;
            }

            // if event list's cycle number is behind, move to the next event
            while(allEvents[channel.first][currentEvent].cycleNumber <
                    macropulseList[currentMacropulse].cycleNumber)
            {
                currentEvent++;
                if(currentEvent>=allEvents[channel.first].size())
                {
                    cout << "Reached end of event list; end event assignment." << endl;
                    endAssignment = true;
                    break;
                }
            }

            if(endAssignment)
            {
                break;
            }

            // move event time to be later than macropulse time
            while((macropulseList[currentMacropulse].macroTime
                        > allEvents[channel.first][currentEvent].completeTime)
                    && (macropulseList[currentMacropulse].cycleNumber==allEvents[channel.first][currentEvent].cycleNumber))
            {
                currentEvent++;
                if(currentEvent>=allEvents[channel.first].size())
                {
                    cout << "Reached end of event list; end event assignment." << endl;
                    endAssignment = true;
                    break;
                }
            }

            if(endAssignment)
            {
                break;
            }

            if(allEvents[channel.first][currentEvent].cycleNumber!=macropulseList[currentMacropulse].cycleNumber)
            {
                // event list has reached the next cycle
                continue;
            }

            // if macropulse's time is more than a full macropulse behind, move to
            // the next macropulse
            if(currentMacropulse >= macropulseList.size()-1)
            {
                break;
            }

            while((macropulseList[currentMacropulse+1].macroTime <
                        allEvents[channel.first][currentEvent].completeTime)
                    && (macropulseList[currentMacropulse+1].cycleNumber==allEvents[channel.first][currentEvent].cycleNumber))
            {
                currentMacropulse++;
                if((int)currentMacropulse>=macropulseList.size())
                {
                    cout << "Reached end of macropulse list; end event assignment." << endl;
                    endAssignment = true;
                    break;
                }

                detectorEvent.eventNo=0;
            }

            if(endAssignment)
            {
                break;
            }

            if(macropulseList[currentMacropulse].cycleNumber!=allEvents[channel.first][currentEvent].cycleNumber)
            {
                // macropulse has reached the next cycle
                currentEvent++;
                continue;
            }

            while((macropulseList[currentMacropulse].macroTime < allEvents[channel.first][currentEvent].completeTime)
                    && ((macropulseList[currentMacropulse+1].macroTime > allEvents[channel.first][currentEvent].completeTime)
                        || (macropulseList[currentMacropulse+1].cycleNumber > allEvents[channel.first][currentEvent].cycleNumber)))
            {
                detectorEvent.cycleNumber = allEvents[channel.first][currentEvent].cycleNumber;
                detectorEvent.completeTime = allEvents[channel.first][currentEvent].completeTime;
                detectorEvent.fineTime = allEvents[channel.first][currentEvent].fineTime;
                detectorEvent.sgQ = allEvents[channel.first][currentEvent].sgQ;
                detectorEvent.lgQ = allEvents[channel.first][currentEvent].lgQ;
                detectorEvent.waveform = allEvents[channel.first][currentEvent].waveform;

                detectorEvent.macroTime = macropulseList[currentMacropulse].macroTime;
                detectorEvent.macroNo = macropulseList[currentMacropulse].macroNo;
                detectorEvent.targetPos = macropulseList[currentMacropulse].targetPos;

                if(channel.second == config.analysis.MONITOR_TREE_NAME)
                {
                    macropulseList[currentMacropulse].numberOfMonitorsInMacro++;
                }

                else if(channel.second == config.cs.DETECTOR_NAMES[0])
                {
                    macropulseList[currentMacropulse].numberOfEventsInMacro++;
                }

                outputTree->Fill();

                detectorEvent.eventNo++;
                currentEvent++;

                if(allEvents[channel.first][currentEvent].cycleNumber!=macropulseList[currentMacropulse].cycleNumber)
                {
                    break;
                }

                if((int)currentEvent>=allEvents[channel.first].size())
                {
                    cout << "Reached end of event list; end event assignment." << endl;
                    endAssignment = true;
                    break;
                }

                if(currentEvent%10000==0)
                {
                    cout << "Assigned " << currentEvent << " events to macropulses...\r";
                    fflush(stdout);
                }
            }

            if(endAssignment)
            {
                break;
            }
        }

        logFile << endl;
        logFile << "For \"" << channel.second << "\" channel:" << endl;
        logFile << "fraction of events successfully assigned to macropulses =  "
            << (double)(outputTree->GetEntries())/allEvents[channel.first].size() << endl;

        outputTree->Write();

    }

    outputFile->Close();
    inputFile->Close();

    return 0;
}
