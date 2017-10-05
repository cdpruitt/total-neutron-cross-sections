/******************************************************************************
  identifyMacropulses.cpp 
 ******************************************************************************/
// This method takes two ROOT Trees as input: one containing all macropulse
// starting time events, and one containing target changer positions at the
// start of each macropulse.  It calculates:
//
// 1) The correct starting time for each macropulse (including a fine time
// correction and, if needed, an extended time corrections).
//
// 2) The target changer position for each macropulse.
//
// The resulting macropulse events are written to another ROOT Tree.

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

#include "../include/identifyMacropulses.h" // declarations of functions used to assign times and macropulses to events

#include "../include/softwareCFD.h"

#include "../include/config.h"

extern Config config;

using namespace std;

// Use the lgQ from the target changer to determine the target position
int assignTargetPos(int lgQ)
{
    for(int i=0; (size_t)i<config.targetConfig.TARGET_GATES.size(); i++)
    {
        if (lgQ>=config.targetConfig.TARGET_GATES[i].first && lgQ<=config.targetConfig.TARGET_GATES[i].second)
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

int identifyMacropulses(
        string inputFileName,
        string inputTreeName,
        string outputFileName,
        string macropulseTreeName)
{
    TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
    if (!inputFile)
    {
        cerr << "Failed to open " << inputFileName << ". Please check that the file exists" << endl;
        return(1);
    }

    TTree* inputTree = (TTree*)inputFile->Get(inputTreeName.c_str());
    if(!inputTree)
    {
        cerr << "Error: couldn't find " << inputTreeName << " tree in "
            << inputFileName << " when attempting to identify macropulses." << endl;
        cerr << "Please check that the tree exists. " << endl;
        return(1);
    }

    // create structs for holding macropulse event data, and link it to the trees
    unsigned int chNo;
    unsigned int cycleNumber;
    double completeTime;
    inputTree->SetBranchAddress("chNo",&chNo);
    inputTree->SetBranchAddress("cycleNumber",&cycleNumber);
    inputTree->SetBranchAddress("completeTime",&completeTime);

    struct MacrotimeEvent
    {
        unsigned int cycleNumber = 0;
        unsigned int macroNo = 0;
        double completeTime = 0;
    } macrotimeEvent;

    struct TargetChangerEvent
    {
        unsigned int cycleNumber = 0;
        double completeTime = 0;
        unsigned int lgQ = 0;
    } targetChangerEvent;

    inputTree->SetBranchAddress("lgQ",&targetChangerEvent.lgQ);

    vector<MacrotimeEvent> macrotimeList;
    vector<TargetChangerEvent> targetChangerList;

    long unsigned int inputTreeEntries = inputTree->GetEntries();
    long unsigned int currentTreeEntry = 0;

    while(currentTreeEntry<inputTreeEntries)
    {
        inputTree->GetEntry(currentTreeEntry);

        if(chNo==0)
        {
            targetChangerEvent.cycleNumber = cycleNumber;
            targetChangerEvent.completeTime = completeTime;
            targetChangerList.push_back(TargetChangerEvent(targetChangerEvent));
        }

        else if(chNo==1)
        {
            macrotimeEvent.cycleNumber = cycleNumber;
            macrotimeEvent.completeTime = completeTime;
            macrotimeList.push_back(MacrotimeEvent(macrotimeEvent));
            macrotimeEvent.macroNo++;
        }

        if(currentTreeEntry%10000==0)
        {
            cout << "Processed " << currentTreeEntry << " events from raw data file looking for macropulse events.\r";
            fflush(stdout);
        }

        currentTreeEntry++;
    }

    cout << endl << "Finished processing events from raw data file. Total macropulses recovered = "
        << macrotimeList.size() << endl;

    // all macrotime and target changer events have been identified;
    // time to combine them into full-fledged macropulse events
    // and write them out to a ROOT tree

    TFile* outputFile = new TFile(outputFileName.c_str(),"CREATE");
    if(!outputFile)
    {
        cerr << "Error: could not create " << outputFileName << " (does it already exist?)." << endl;
        return(1);
    }

    struct MacropulseEvent
    {
        unsigned int cycleNumber = 0;
        unsigned int macroNo = 0;
        double macroTime = 0;
        unsigned int targetPos = 0;
        unsigned int lgQ = 0;
    } macropulseEvent;

    TTree* macropulseTree = new TTree(macropulseTreeName.c_str(),"");

    macropulseTree->Branch("macroTime",&macropulseEvent.macroTime,"macroTime/d");
    macropulseTree->Branch("cycleNumber",&macropulseEvent.cycleNumber,"cycleNumber/i");
    macropulseTree->Branch("macroNo",&macropulseEvent.macroNo,"macroNo/i");
    macropulseTree->Branch("lgQ",&macropulseEvent.lgQ,"lgQ/i");
    macropulseTree->Branch("targetPos",&macropulseEvent.targetPos,"targetPos/i");

    unsigned int currentMacrotimeEntry = 0;
    unsigned int currentTargetChangerEntry = 0;

    bool endAssignment = false;

    while(currentMacrotimeEntry<macrotimeList.size())
    {
        // if macrotimeList's cycle number is behind, move to the next macrotime
        // event
        while((macrotimeList[currentMacrotimeEntry].cycleNumber <
               targetChangerList[currentTargetChangerEntry].cycleNumber)
             )
        {
            currentMacrotimeEntry++;
            if(currentMacrotimeEntry>=macrotimeList.size())
            {
                cout << "Reached end of macrotime list; end macropulse identification." << endl;
                endAssignment = true;
                break;
            }
        }

        // if targetChangerList's cycle number is behind, move to the next
        // target changer event
        while((macrotimeList[currentMacrotimeEntry].cycleNumber >
               targetChangerList[currentTargetChangerEntry].cycleNumber)
             )
        {
            currentTargetChangerEntry++;
            if(currentTargetChangerEntry>=targetChangerList.size())
            {
                cout << "Reached end of target changer list; end macropulse identification." << endl;
                endAssignment = true;
                break;
            }
        }

        if(endAssignment)
        {
            break;
        }

        // macrotime and target changer events are in the same cycle
        // if macrotimeList's time is behind, move to the next macrotime
        // event
        while((macrotimeList[currentMacrotimeEntry].completeTime+20 <
                 targetChangerList[currentTargetChangerEntry].completeTime)
             )
        {
            currentMacrotimeEntry++;
            if(currentMacrotimeEntry>=macrotimeList.size())
            {
                cout << "Reached end of macrotime list; end macropulse identification." << endl;
                endAssignment = true;
                break;
            }
        }

        // if targetChangerList's time is behind, move to the next target
        // changer event
        while((macrotimeList[currentMacrotimeEntry].completeTime >
               targetChangerList[currentTargetChangerEntry].completeTime+20)
             )
        {
            currentTargetChangerEntry++;
            if(currentTargetChangerEntry>=targetChangerList.size())
            {
                cout << "Reached end of target changer list; end macropulse identification." << endl;
                endAssignment = true;
                break;
            }
        }

        if((macrotimeList[currentMacrotimeEntry].cycleNumber == 
                    targetChangerList[currentTargetChangerEntry].cycleNumber)
                && (abs(macrotimeList[currentMacrotimeEntry].completeTime - 
                        targetChangerList[currentTargetChangerEntry].completeTime)<20)
          )
        {
            // found a match between the target changer and the macrotime lists
            macropulseEvent.lgQ = targetChangerList[currentTargetChangerEntry].lgQ;
            macropulseEvent.targetPos = assignTargetPos(macropulseEvent.lgQ);
        }

        else
        {
            // failed to find a match; assign a target position of 0
            macropulseEvent.lgQ = 1;
            macropulseEvent.targetPos = 0;
        }

        macropulseEvent.cycleNumber = macrotimeList[currentMacrotimeEntry].cycleNumber;
        macropulseEvent.macroNo = macrotimeList[currentMacrotimeEntry].macroNo;
        macropulseEvent.macroTime = macrotimeList[currentMacrotimeEntry].completeTime;

        macropulseTree->Fill();

        if(currentMacrotimeEntry%100==0)
        {
            cout << "Identified " << currentMacrotimeEntry << " macropulses...\r";
            fflush(stdout);
        }

        currentMacrotimeEntry++;
    }

    // clean up
    macropulseTree->Write();
    outputFile->Close();

    inputFile->Close();

    return 0;
}
