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
    for(int i=0; (size_t)i<config.target.TARGET_GATES.size(); i++)
    {
        if (lgQ>=config.target.TARGET_GATES[i].first && lgQ<=config.target.TARGET_GATES[i].second)
        {
            // lgQ fits within this gate
            return i; // target positions start from 1
        }
    }

    return -1;
}

int identifyMacropulses(
        string inputFileName,
        string outputFileName,
        ofstream& logFile)
{
    // check to see if output file already exists; if so, exit
    ifstream f(outputFileName);

    if(f.good())
    {
         cout << outputFileName << " already exists; skipping event assignment to macropulses." << endl;
         logFile << outputFileName << " already exists; skipping event assignment to macropulses." << endl;
         return 2;
    }

    f.close();

    TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
    if (!inputFile)
    {
        cerr << "Failed to open " << inputFileName << ". Please check that the file exists" << endl;
        return 1;
    }

    TTree* inputTree = (TTree*)inputFile->Get(config.analysis.DPP_TREE_NAME.c_str());
    if(!inputTree)
    {
        cerr << "Error: couldn't find " << config.analysis.DPP_TREE_NAME << " tree in "
            << inputFileName << " when attempting to identify macropulses." << endl;
        cerr << "Please check that the tree exists. " << endl;
        return 1;
    }

    // create structs for holding macropulse event data, and link it to the trees

    struct MacrotimeEvent
    {
        unsigned int cycleNumber = 0;
        unsigned int macroNo = 0;
        double completeTime = 0;
        vector<int> waveform;
    } macrotimeEvent;

    struct TargetChangerEvent
    {
        unsigned int cycleNumber = 0;
        double completeTime = 0;
        unsigned int lgQ = 0;
    } targetChangerEvent;

    unsigned int chNo;
    unsigned int cycleNumber;
    double completeTime;
    vector<int>* waveformPointer = 0;

    inputTree->SetBranchAddress("chNo",&chNo);
    inputTree->SetBranchAddress("cycleNumber",&cycleNumber);
    inputTree->SetBranchAddress("completeTime",&completeTime);
    inputTree->SetBranchAddress("lgQ",&targetChangerEvent.lgQ);
    inputTree->SetBranchAddress("waveform",&waveformPointer);

    vector<MacrotimeEvent> macroTimeList;
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
            macrotimeEvent.waveform = *waveformPointer;
            macroTimeList.push_back(MacrotimeEvent(macrotimeEvent));
            macrotimeEvent.macroNo++;
        }

        if(currentTreeEntry%10000==0)
        {
            cout << "Processed " << currentTreeEntry << " events from raw data file looking for macropulse events.\r";
            fflush(stdout);
        }

        currentTreeEntry++;
    }

    cout << "Finished processing events from raw data file. Total macropulses recovered = "
        << macroTimeList.size() << endl;

    logFile << endl << "Finished processing events from raw data file. Total macropulses recovered = "
        << macroTimeList.size() << endl;

    // all macrotime and target changer events have been identified;
    // time to combine them into full-fledged macropulse events
    // and write them out to a ROOT tree
    TFile* outputFile = new TFile(outputFileName.c_str(),"CREATE");
    if(!outputFile)
    {
        cerr << "Error: could not create " << outputFileName << " (does it already exist?)." << endl;
        return 1;
    }

    MacropulseEvent macropulseEvent;
    waveformPointer = new vector<int>;

    TTree* macropulseTree = new TTree(config.analysis.MACROPULSE_TREE_NAME.c_str(),"");

    macropulseTree->Branch("macroTime",&macropulseEvent.macroTime,"macroTime/d");
    macropulseTree->Branch("cycleNumber",&macropulseEvent.cycleNumber,"cycleNumber/i");
    macropulseTree->Branch("macroNo",&macropulseEvent.macroNo,"macroNo/i");
    macropulseTree->Branch("lgQ",&macropulseEvent.lgQ,"lgQ/i");
    macropulseTree->Branch("targetPos",&macropulseEvent.targetPos,"targetPos/I");
    macropulseTree->Branch("waveform",&(macropulseEvent.waveform));

    unsigned int currentMacrotimeEntry = 0;
    unsigned int currentTargetChangerEntry = 0;

    bool endAssignment = false;

    while(currentMacrotimeEntry<macroTimeList.size())
    {
        // if macroTimeList's cycle number is behind, move to the next macrotime
        // event
        while((macroTimeList[currentMacrotimeEntry].cycleNumber <
               targetChangerList[currentTargetChangerEntry].cycleNumber)
             )
        {
            currentMacrotimeEntry++;
            if(currentMacrotimeEntry>=macroTimeList.size())
            {
                cout << "Reached end of macrotime list; end macropulse identification." << endl;
                endAssignment = true;
                break;
            }
        }

        // if targetChangerList's cycle number is behind, move to the next
        // target changer event
        while((macroTimeList[currentMacrotimeEntry].cycleNumber >
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
        // if macroTimeList's time is behind, move to the next macrotime
        // event
        while((macroTimeList[currentMacrotimeEntry].completeTime+20 <
                 targetChangerList[currentTargetChangerEntry].completeTime)
             )
        {
            currentMacrotimeEntry++;
            if(currentMacrotimeEntry>=macroTimeList.size())
            {
                cout << "Reached end of macrotime list; end macropulse identification." << endl;
                endAssignment = true;
                break;
            }
        }

        // if targetChangerList's time is behind, move to the next target
        // changer event
        while((macroTimeList[currentMacrotimeEntry].completeTime >
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

        if((macroTimeList[currentMacrotimeEntry].cycleNumber == 
                    targetChangerList[currentTargetChangerEntry].cycleNumber)
                && (abs(macroTimeList[currentMacrotimeEntry].completeTime - 
                        targetChangerList[currentTargetChangerEntry].completeTime)<20)
          )
        {
            // found a match between the target changer and the macrotime lists
            macropulseEvent.lgQ = targetChangerList[currentTargetChangerEntry].lgQ;
            macropulseEvent.targetPos = assignTargetPos(macropulseEvent.lgQ);

            if(macropulseEvent.targetPos < 0)
            {
                logFile << "Target changer assignment error at macropulse "
                    << macroTimeList[currentMacrotimeEntry].macroNo
                    << ": lgQ was " << macropulseEvent.lgQ << endl;
                currentMacrotimeEntry++;
                continue;
            }
        }

        else
        {
            // failed to find a match; discard this macropulse
            currentMacrotimeEntry++;
            continue;
        }

        macropulseEvent.cycleNumber = macroTimeList[currentMacrotimeEntry].cycleNumber;
        macropulseEvent.macroNo = macroTimeList[currentMacrotimeEntry].macroNo;
        macropulseEvent.macroTime = macroTimeList[currentMacrotimeEntry].completeTime;
        macropulseEvent.waveform = macroTimeList[currentMacrotimeEntry].waveform;

        macropulseTree->Fill();

        if(currentMacrotimeEntry%100==0)
        {
            cout << "Identified " << currentMacrotimeEntry << " macropulses...\r";
            fflush(stdout);
        }

        currentMacrotimeEntry++;
    }

    cout << "Finished macropulse identification. Total number of macropulses identified = "
        << macropulseTree->GetEntries() << endl;

    logFile << endl;
    logFile << "Total number of macropulse times identified: " << macroTimeList.size() << "." << endl;
    logFile << "Total number of target changer events identified: " << targetChangerList.size() << "." << endl;

    logFile << "Fraction of macropulse times mated to a target changer event: "
        << macropulseTree->GetEntries() << endl;

    // clean up
    macropulseTree->Write();
    outputFile->Close();

    inputFile->Close();

    return 0;
}
