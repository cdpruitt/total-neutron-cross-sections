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

int assignEventsToMacropulses(string inputFileName, string outputFileName, ofstream& logFile, pair<unsigned int, string> channel)
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

    /**************************************************************************/
    // open macropulse tree
    TFile* outputFile = new TFile(outputFileName.c_str(), "UPDATE");
    if(!outputFile)
    {
        cerr << "Failed to open " << outputFileName << ". Please check that the file exists" << endl;
        exit(1);
    }

    TTree* macropulseTree = (TTree*)outputFile->Get(config.analysis.MACROPULSE_TREE_NAME.c_str());
    if(!macropulseTree)
    {
        cerr << "Error: couldn't find " << config.analysis.MACROPULSE_TREE_NAME << " tree in "
             << inputFileName << " when attempting to assign macropulses." << endl;
        cerr << "Please check that the tree exists. " << endl;
        return(1);
    }
    /**************************************************************************/
 
    /**************************************************************************/
    // create struct for holding event data and link it to trees
    DetectorEvent detectorEvent;
    
    unsigned int chNo;
    vector<int>* waveformPointer = 0;

    inputTree->SetBranchAddress("chNo",&chNo);
    inputTree->SetBranchAddress("cycleNumber",&detectorEvent.cycleNumber);
    inputTree->SetBranchAddress("completeTime",&detectorEvent.completeTime);
    inputTree->SetBranchAddress("fineTime",&detectorEvent.fineTime);
    inputTree->SetBranchAddress("sgQ",&detectorEvent.sgQ);
    inputTree->SetBranchAddress("lgQ",&detectorEvent.lgQ);
    inputTree->SetBranchAddress("waveform",&waveformPointer);

    /**************************************************************************/

    vector<DetectorEvent> eventList;

    long unsigned int inputTreeEntries = inputTree->GetEntries();
    long unsigned int currentTreeEntry = 0;

    while(currentTreeEntry<inputTreeEntries)
    {
        inputTree->GetEntry(currentTreeEntry);

        if(chNo==channel.first)
        {
            DetectorEvent de = DetectorEvent(detectorEvent);
            de.waveform = *waveformPointer;
            eventList.push_back(de);
            detectorEvent.eventNo++;
        }

        if(currentTreeEntry%10000==0)
        {
            cout << "Read " << currentTreeEntry << " \"" << config.analysis.DPP_TREE_NAME << "\" events...\r";
            fflush(stdout);
        }

        currentTreeEntry++;
    }

    cout << "Finished reading \"" << channel.second << "\" events from raw data file. Total events recovered = "
        << eventList.size() << endl;

    /**************************************************************************/
    // create a new tree for holding sorted detector events
    TTree* outputTree = new TTree(channel.second.c_str(),"");
    if(!outputTree)
    {
        cerr << "Error: couldn't create detector tree " << channel.second << " when attempting to assign events to macropulses. Exiting... " << endl;
        return 1;
    }
    /**************************************************************************/

    unsigned int macropulseCycleNumber;
    unsigned int macroNo;
    double macroTime;
    int targetPos;

    macropulseTree->SetBranchAddress("cycleNumber",&macropulseCycleNumber);
    macropulseTree->SetBranchAddress("macroNo",&macroNo);
    macropulseTree->SetBranchAddress("macroTime",&macroTime);
    macropulseTree->SetBranchAddress("targetPos",&targetPos);

    outputTree->Branch("macroTime",&detectorEvent.macroTime,"macroTime/D");
    outputTree->Branch("completeTime",&detectorEvent.completeTime,"completeTime/d");
    outputTree->Branch("fineTime",&detectorEvent.fineTime,"fineTime/d");

    outputTree->Branch("macroNo",&detectorEvent.macroNo,"macroNo/i");
    outputTree->Branch("cycleNumber",&detectorEvent.cycleNumber,"cycleNumber/i");
    outputTree->Branch("targetPos",&detectorEvent.targetPos,"targetPos/I");
    outputTree->Branch("eventNo",&detectorEvent.eventNo,"eventNo/i");
    outputTree->Branch("sgQ",&detectorEvent.sgQ,"sgQ/i");
    outputTree->Branch("lgQ",&detectorEvent.lgQ,"lgQ/i");
    outputTree->Branch("waveform",&detectorEvent.waveform);

    unsigned int currentEvent = 0;
    unsigned int currentMacropulseEntry = 0;
    macropulseTree->GetEntry(currentMacropulseEntry);

    bool endAssignment = false;

    unsigned int currentCycleNumber = 0;
    unsigned int maxCycleNumber = eventList.back().cycleNumber;

    const double MACRO_PERIOD = pow(10,6)*8; // in ns

    while(currentCycleNumber<=maxCycleNumber)
    {
        // if the macropulse's cycle number is behind, move to the next macropulse
        while(macropulseCycleNumber <
                currentCycleNumber)
        {
            currentMacropulseEntry++;
            if(currentMacropulseEntry>=macropulseTree->GetEntries())
            {
                cout << "Reached end of macrotime list; end event assignment." << endl;
                endAssignment = true;
                break;
            }

            macropulseTree->GetEntry(currentMacropulseEntry);
            detectorEvent.eventNo=0;
        }

        // if event list's cycle number is behind, move to the next event
        while(eventList[currentEvent].cycleNumber <
                currentCycleNumber)
        {
            currentEvent++;
            if(currentEvent>=eventList.size())
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

        // Both the macropulse and the current event are in the current cycle

        // if macropulse's time is more than a full macropulse behind, move to
        // the next macropulse
        while(macroTime + MACRO_PERIOD <
                 eventList[currentEvent].completeTime
           && macropulseCycleNumber==currentCycleNumber)
        {
            currentMacropulseEntry++;
            if(currentMacropulseEntry>=macropulseTree->GetEntries())
            {
                cout << "Reached end of macropulse list; end event assignment." << endl;
                endAssignment = true;
                break;
            }

            macropulseTree->GetEntry(currentMacropulseEntry);
            detectorEvent.eventNo=0;
        }

        if(endAssignment)
        {
            break;
        }

        if(macropulseCycleNumber!=currentCycleNumber)
        {
            // macropulse has reached the next cycle
            currentCycleNumber++;
            continue;
        }

        while((macroTime > eventList[currentEvent].completeTime)
           && (currentCycleNumber==eventList[currentEvent].cycleNumber))
        {
            currentEvent++;
            if(currentEvent>eventList.size())
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

        if(eventList[currentEvent].cycleNumber!=currentCycleNumber)
        {
            // event list has reached the next cycle
            currentCycleNumber++;
            continue;
        }

        while((macroTime < eventList[currentEvent].completeTime)
           && (macroTime + MACRO_PERIOD > eventList[currentEvent].completeTime)
           && (eventList[currentEvent].cycleNumber==currentCycleNumber))
        {
            detectorEvent.cycleNumber = eventList[currentEvent].cycleNumber;
            detectorEvent.completeTime = eventList[currentEvent].completeTime;
            detectorEvent.fineTime = eventList[currentEvent].fineTime;
            detectorEvent.sgQ = eventList[currentEvent].sgQ;
            detectorEvent.lgQ = eventList[currentEvent].lgQ;
            detectorEvent.waveform = eventList[currentEvent].waveform;

            detectorEvent.macroTime = macroTime;
            detectorEvent.macroNo = macroNo;
            detectorEvent.targetPos = targetPos;

            outputTree->Fill();
            detectorEvent.eventNo++;

            currentEvent++;
            if(currentEvent>eventList.size())
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

        if(eventList[currentEvent].cycleNumber!=currentCycleNumber)
        {
            // event list has reached the next cycle
            currentCycleNumber++;
            continue;
        }
    }

    logFile << endl;
    logFile << "For \"" << channel.second << "\" channel:" << endl;
    logFile << "fraction of events successfully assigned to macropulses =  "
        << (double)(outputTree->GetEntries())/eventList.size() << endl;

    outputTree->Write();
    outputFile->Close();
    inputFile->Close();

    return 0;
}
