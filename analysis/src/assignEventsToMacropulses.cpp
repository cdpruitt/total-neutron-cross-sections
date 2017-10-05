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

extern Config config;

using namespace std;

// create struct for holding detector event data and link to trees
struct DetectorEvent
{
    double macroTime = 0;
    unsigned int macroNo = 0;
    unsigned int targetPos = 0;

    unsigned int cycleNumber = 0;
    double completeTime = 0;
    unsigned int timetag = 0;
    unsigned int extTime = 0;
    double fineTime = 0;
    unsigned int eventNo = 0;
    unsigned int sgQ = 0;
    unsigned int lgQ = 0;
    unsigned int baseline = 0;
    vector<int>* waveform = new vector<int>;
};

int assignEventsToMacropulses(string inputFileName, string inputTreeName, string outputFileName, string macropulseTreeName, unsigned int channelNo, string outputTreeName)
{
    /**************************************************************************/
    // open input detector tree
    TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
    if(!inputFile)
    {
        cerr << "Failed to open " << inputFileName << ". Please check that the file exists" << endl;
        exit(1);
    }
    cout << inputFileName << " opened successfully. Start reading events..." << endl;

    TTree* inputTree = (TTree*)inputFile->Get(inputTreeName.c_str());
    if(!inputTree)
    {
        cerr << "Error: couldn't find detector tree " << inputTreeName << " in "
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

    TTree* macropulseTree = (TTree*)outputFile->Get(macropulseTreeName.c_str());
    if(!macropulseTree)
    {
        cerr << "Error: couldn't find " << macropulseTreeName << " tree in "
             << inputFileName << " when attempting to assign macropulses." << endl;
        cerr << "Please check that the tree exists. " << endl;
        return(1);
    }
    cout << macropulseTreeName << " opened successfully." << endl;
    /**************************************************************************/
 
    /**************************************************************************/
    // create struct for holding event data and link it to trees
    DetectorEvent detectorEvent;
    
    unsigned int chNo;
    inputTree->SetBranchAddress("chNo",&chNo);
    inputTree->SetBranchAddress("cycleNumber",&detectorEvent.cycleNumber);
    inputTree->SetBranchAddress("completeTime",&detectorEvent.completeTime);
    inputTree->SetBranchAddress("fineTime",&detectorEvent.fineTime);
    inputTree->SetBranchAddress("sgQ",&detectorEvent.sgQ);
    inputTree->SetBranchAddress("lgQ",&detectorEvent.lgQ);
    inputTree->SetBranchAddress("waveform",&detectorEvent.waveform);

    /**************************************************************************/

    vector<DetectorEvent> eventList;

    long unsigned int inputTreeEntries = inputTree->GetEntries();
    long unsigned int currentTreeEntry = 0;

    while(currentTreeEntry<inputTreeEntries)
    {
        inputTree->GetEntry(currentTreeEntry);

        if(chNo==channelNo)
        {
            eventList.push_back(DetectorEvent(detectorEvent));
            detectorEvent.eventNo++;
        }

        if(currentTreeEntry%10000==0)
        {
            cout << "Processed " << currentTreeEntry << " events from raw data file while idenfying detector events.\r";
            fflush(stdout);
        }

        currentTreeEntry++;
    }

    cout << endl << "Finished processing events from raw data file. Total events recovered = "
        << eventList.size() << endl;

    /**************************************************************************/
    // create a new tree for holding sorted detector events
    TTree* outputTree = new TTree(outputTreeName.c_str(),"");
    if(!outputTree)
    {
        cerr << "Error: couldn't create detector tree " << outputTreeName << " when attempting to assign events to macropulses. Exiting... " << endl;
        return 1;
    }
    /**************************************************************************/

    unsigned int macropulseCycleNumber;
    unsigned int macroNo;
    double macroTime;
    unsigned int targetPos;

    macropulseTree->SetBranchAddress("cycleNumber",&macropulseCycleNumber);
    macropulseTree->SetBranchAddress("macroNo",&macroNo);
    macropulseTree->SetBranchAddress("macroTime",&macroTime);
    macropulseTree->SetBranchAddress("targetPos",&targetPos);

    outputTree->Branch("macroTime",&detectorEvent.macroTime,"macroTime/d");
    outputTree->Branch("macroNo",&detectorEvent.macroNo,"macroNo/i");
    outputTree->Branch("cycleNumber",&detectorEvent.cycleNumber,"cycleNumber/i");
    outputTree->Branch("targetPos",&detectorEvent.targetPos,"targetPos/i");

    outputTree->Branch("completeTime",&detectorEvent.completeTime,"completeTime/d");
    outputTree->Branch("fineTime",&detectorEvent.fineTime,"fineTime/d");
    outputTree->Branch("eventNo",&detectorEvent.eventNo,"eventNo/i");
    outputTree->Branch("sgQ",&detectorEvent.sgQ,"sgQ/i");
    outputTree->Branch("lgQ",&detectorEvent.lgQ,"lgQ/i");
    outputTree->Branch("waveform",&detectorEvent.waveform);

    unsigned int currentEvent = 0;
    unsigned int currentMacropulseEntry = 0;

    macropulseTree->GetEntry(currentMacropulseEntry);

    bool endAssignment = false;

    cout << "Assigning " << outputTreeName << " events to macropulses..." << endl;
    while(currentEvent<eventList.size())
    {
        // if the macropulse's cycle number is behind, move to the next macropulse
        while(macropulseCycleNumber <
                eventList[currentEvent].cycleNumber)
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
        while(macropulseCycleNumber >
               eventList[currentEvent].cycleNumber)
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

        // if detector event's time is behind, move to the next detector event
        while((macroTime >
                    eventList[currentEvent].completeTime)
             )
        {
            currentEvent++;
            if(currentEvent>=eventList.size())
            {
                cout << "Reached end of event list; end event assignement." << endl;
                endAssignment = true;
                break;
            }
        }

        // macropulse and event are in the same cycle
        // if macropulse's time is behind, move to the next macropulse
        // event
        while((macroTime + config.facilityConfig.MACRO_LENGTH <
                 eventList[currentEvent].completeTime)
             )
        {
            currentMacropulseEntry++;
            if(currentMacropulseEntry>=macropulseTree->GetEntries())
            {
                cout << "Reached end of macropulse list; end event assignement." << endl;
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

        if(!(macropulseCycleNumber == 
                eventList[currentEvent].cycleNumber
            )
          )
        {
            // failed to find the correct macro for this event; continue to next detector event
            continue;
        }

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

        if(currentEvent%10000==0)
        {
            cout << "Assigned " << currentEvent << " events to macropulses...\r";
            fflush(stdout);
        }

        currentEvent++;
    }

    outputTree->Write();
    outputFile->Close();
    inputFile->Close();

    return 0;
}
