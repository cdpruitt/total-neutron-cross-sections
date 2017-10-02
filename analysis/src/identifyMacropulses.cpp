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

int identifyMacropulses(
        string inputFileName,
        string macropulseTimeTreeName,
        string targetPositionTreeName,
        string outputFileName,
        string outputTreeName)
{
    TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
    if (!inputFile)
    {
        cerr << "Failed to open " << inputFileName << ". Please check that the file exists" << endl;
        return(1);
    }

    TTree* macropulseTimeTree = (TTree*)inputFile->Get(macropulseTimeTreeName.c_str());
    if(!macropulseTimeTree)
    {
        cerr << "Error: couldn't find " << macropulseTimeTreeName << " tree in "
             << inputFileName << " when attempting to identify macropulses." << endl;
        cerr << "Please check that the tree exists. " << endl;
        return(1);
    }

    TTree* targetPositionTree = (TTree*)inputFile->Get(targetPositionTreeName.c_str());
    if(!targetPositionTree)
    {
        cerr << "Error: couldn't find " << targetPositionTreeName << " tree in "
             << inputFileName << " when attempting to identify macropulses." << endl;
        cerr << "Please check that the tree exists. " << endl;
        return(1);
    }

    TFile* outputFile = new TFile(outputFileName.c_str(),"CREATE");

    if(!outputFile)
    {
        cerr << "Error: could not create " << outputFileName << " (does it already exist?)." << endl;
        return(1);
    }

    TTree* outputTree = new TTree(outputTreeName.c_str(),"");

    // create struct for holding macropulse event data, and link it to the trees
    struct MacropulseEvent
    {
        unsigned int timeTimetag = 0;
        unsigned int timeExtTime = 0;
        vector<int>* timeWaveform = new vector<int>;

        unsigned int targetTimetag = 0;
        unsigned int targetExtTime = 0;
        unsigned int lgQ = 0;
        vector<int>* targetWaveform = new vector<int>;

        double macroTime = 0;
        double fineTime = 0;
        unsigned int macroNo = 0;
        unsigned int targetPos = 0;
        unsigned int modeChange = 0;

    } macropulseEvent;

    macropulseTimeTree->SetBranchAddress("timetag",&macropulseEvent.timeTimetag);
    macropulseTimeTree->SetBranchAddress("extTime",&macropulseEvent.timeExtTime);
    macropulseTimeTree->SetBranchAddress("waveform",&macropulseEvent.timeWaveform);

    targetPositionTree->SetBranchAddress("timetag",&macropulseEvent.targetTimetag);
    targetPositionTree->SetBranchAddress("extTime",&macropulseEvent.targetExtTime);
    targetPositionTree->SetBranchAddress("lgQ",&macropulseEvent.lgQ);
    targetPositionTree->SetBranchAddress("waveform",&macropulseEvent.targetWaveform);

    outputTree->Branch("macroTime",&macropulseEvent.macroTime,"macroTime/D");
    outputTree->Branch("macroNo",&macropulseEvent.macroNo,"macroNo/i");
    outputTree->Branch("modeChange",&macropulseEvent.modeChange,"modeChange/i");
    outputTree->Branch("targetPos",&macropulseEvent.targetPos,"targetPos/i");
    outputTree->Branch("lgQ",&macropulseEvent.lgQ,"lgQ/i");
    outputTree->Branch("timeWaveform",&macropulseEvent.timeWaveform);
    outputTree->Branch("targetWaveform",&macropulseEvent.targetWaveform);

    // loop through all events in the input trees to identify good macropulses
    unsigned int prevExtTime = 0;
    unsigned int prevMacroTime = 0;

    long unsigned int macropulseTimeTreeEntries = macropulseTimeTree->GetEntries();
    long unsigned int targetPositionTreeEntries = targetPositionTree->GetEntries();

    long unsigned int i = 0;
    long unsigned int entryNumberOffset = 0;

    while(i<targetPositionTreeEntries)
    {
        targetPositionTree->GetEntry(i);
        macropulseTimeTree->GetEntry(i+entryNumberOffset);

        // determine time of target position measurement event
        double targetPositionTime =
            (double)(macropulseEvent.targetExtTime)*pow(2,31)
          + (double)macropulseEvent.targetTimetag;
        targetPositionTime *= experimentalConfig.timeConfig.SAMPLE_PERIOD;

        // determine macropulse start time
        macropulseEvent.macroTime =
            (double)(macropulseEvent.timeExtTime)*pow(2,31)
          + (double)macropulseEvent.timeTimetag;
        macropulseEvent.macroTime *= experimentalConfig.timeConfig.SAMPLE_PERIOD;

        // if macropulse start time and target position measurement time are
        // aligned, assign the target position measurement to the macropulse
        double timeDifferenceBetweenEvents =
            macropulseEvent.macroTime -
            targetPositionTime;

        while(timeDifferenceBetweenEvents < experimentalConfig.timeConfig.MACROPULSE_TARGET_TIME_DIFFERENCE - 15
           || timeDifferenceBetweenEvents > experimentalConfig.timeConfig.MACROPULSE_TARGET_TIME_DIFFERENCE + 15)
        {
            cerr << "Timing mismatch between time of macropulse start and time of target position measurement. Skipping to next macropulse..." << endl;

            cerr << "macropulseNo = " << macropulseEvent.macroNo << endl;
            cerr << "timeDifference = " << timeDifferenceBetweenEvents << endl;

            // move the macropulse time tree forward one entry
            entryNumberOffset++;

            if(i+entryNumberOffset > macropulseTimeTreeEntries)
            {
                cerr << "Reached the end of the macropulse time tree. Ending macropulse identification." << endl;

                outputTree->Write();
                outputFile->Close();

                inputFile->Close();

                return 0;
            }

            macropulseTimeTree->GetEntry(i+entryNumberOffset);

            // determine macropulse start time
            macropulseEvent.macroTime =
                (double)(macropulseEvent.timeExtTime)*pow(2,31)
                + (double)macropulseEvent.timeTimetag;
            macropulseEvent.macroTime *= experimentalConfig.timeConfig.SAMPLE_PERIOD;

            timeDifferenceBetweenEvents =
                macropulseEvent.macroTime -
                targetPositionTime;
        }

        // Discard events with digitizer rollover error (i.e., incrementing extTime before clearing the timetag, near the rollover period at 2^31 bits)

        if(macropulseEvent.timeExtTime > prevExtTime
                && macropulseEvent.timeTimetag > pow(2,31)-1000)
        {
            cerr << "Error: digitizer rollover error. Skipping event. macroNo = " << macropulseEvent.macroNo << ", extTime = " << macropulseEvent.timeExtTime << ", timetag = " << macropulseEvent.timeTimetag << endl;
            cerr << "Skipping to next target changer event..." << endl;

            i++;
            continue;
        }

        macropulseEvent.targetPos = assignTargetPos(macropulseEvent.lgQ);

        // if this is the first event in a DPP mode period, note this
        if(prevMacroTime > macropulseEvent.macroTime)
        {
            macropulseEvent.modeChange = 1;
        }

        else
        {
            macropulseEvent.modeChange = 0;
        }

        outputTree->Fill();

        // increment macropulse number
        macropulseEvent.macroNo++;

        prevMacroTime = macropulseEvent.macroTime;
        prevExtTime = macropulseEvent.timeExtTime;

        if(i%1000==0)
        {
            cout << "processed " << i << " events in target changer tree.\r";
            fflush(stdout);
        }

        i++;
    }

    // clean up
    outputTree->Write();
    outputFile->Close();

    inputFile->Close();

    return 0;
}
