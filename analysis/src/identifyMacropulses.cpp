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
            << "Please check that the tree exists. " << endl;
        return(1);
    }

    TTree* targetPositionTree = (TTree*)inputFile->Get(targetPositionTreeName.c_str());
    if(!targetPositionTree)
    {
        cerr << "Error: couldn't find " << targetPositionTreeName << " tree in "
            << inputFileName << " when attempting to identify macropulses." << endl;
            << "Please check that the tree exists. " << endl;
        return(1);
    }

    MacropulseEvent macropulseEvent;

    macropulseTimeTree->SetBranchAddress("timetag",&macropulseEvent.timeTimetag);
    macropulseTimeTree->SetBranchAddress("extTime",&macropulseEvent.timeExtTime);
    macropulseTimeTree->SetBranchAddress("fineTime",&macropulseEvent.fineTime);
    macropulseTimeTree->SetBranchAddress("waveform",&macropulseEvent.timeWaveform);

    targetPositionTree->SetBranchAddress("timetag",&macropulseEvent.targetTimetag);
    targetPositionTree->SetBranchAddress("extTime",&macropulseEvent.targetExtTime);
    targetPositionTree->SetBranchAddress("lgQ",&macropulseEvent.lgQ);
    targetPositionTree->SetBranchAddress("waveform",&macropulseEvent.targetWaveform);


    TFile* outputFile = new TFile(outputFileName.c_str(),"CREATE");

    if(!outputFile)
    {
        cerr << "Error: could not create " << outputFileName << " (does it already exist?)." << endl;
        return(1);
    }

    TTree* outputTree = new TTree(outputTreeName.c_str(),"");
    tree->Branch("macroTime",&macropulseEvent.macroTime,"macroTime/D");
    tree->Branch("macroNo",&macropulseEvent.macroNo,"macroNo/i");
    tree->Branch("extTime",&macropulseEvent.extTime,"extTime/d");
    tree->Branch("fineTime",&macropulseEvent.fineTime,"fineTime/d");
    tree->Branch("modeChange",&macropulseEvent.modeChange,"modeChange/i");
    tree->Branch("targetPos",&macropulseEvent.targetPos,"targetPos/i");
    tree->Branch("lgQ",&macropulseEvent.lgQ,"lgQ/i");
    tree->Branch("timeWaveform",&macropulseEvent.timeWaveform);
    tree->Branch("targetWaveform",&macropulseEvent.targetWaveform);

    double prevTimetag = 0;

    long totalEntries = macropulseTimeTree->GetEntries();

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

    TargetChangerEvent tcEvent;

    branchTargetChanger(outputTree, tcEvent);

    for(int j=0; j<totalEntries; j++)
    {
        macropulseTimeTree->GetEntry(j);

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

        outputTree->Fill();

        // update macropulse counters
        tcEvent.macroNo++;

        prevTimetag = tcEvent.macroTime;

        if(j%1000==0)
        {
            cout << "processed " << j << " events in target changer tree.\r";
            fflush(stdout);
        }
    }

    outputTree->Write();
    outputFile->Close();
}
