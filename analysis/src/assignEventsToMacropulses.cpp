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
#include "../include/softwareCFD.h" // to calculate fine times of events
#include "../include/experimentalConfig.h"

extern ExperimentalConfig experimentalConfig;

using namespace std;

// create struct for holding detector event data and link to trees
struct DetectorEvent
{
    double macroTime = 0;
    unsigned int macroNo = 0;
    unsigned int targetPos = 0;

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

void determineDetectorEventTime(DetectorEvent& detectorEvent)
{
    detectorEvent.completeTime =
            (double)(detectorEvent.extTime)*pow(2,31)
            + (double)detectorEvent.timetag;

    // add detector fine time, if it can be calculated
        detectorEvent.fineTime = calculateCFDTime(
                detectorEvent.waveform,
                detectorEvent.baseline,
                experimentalConfig.timeConfig.CFD_FRACTION,
                experimentalConfig.timeConfig.CFD_DELAY);

        if(detectorEvent.fineTime>=0)
            {
                // recovered a good fine time for this event
                detectorEvent.completeTime += detectorEvent.fineTime;
                detectorEvent.completeTime -= experimentalConfig.timeConfig.FINE_TIME_OFFSET;
            }

        // change units to ns from samples
            detectorEvent.completeTime *= experimentalConfig.timeConfig.SAMPLE_PERIOD;

            // correct for the cable delay of the macropulse timing signal
            detectorEvent.completeTime += experimentalConfig.timeConfig.SUMMED_DETECTOR_TIME_OFFSET;
}

int assignEventsToMacropulses(string inputFileName, string detectorTreeName, string outputFileName, string macropulseTreeName)
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

    TTree* inputDetectorTree = (TTree*)inputFile->Get(detectorTreeName.c_str());
    if(!inputDetectorTree)
    {
        cerr << "Error: couldn't find detector tree " << detectorTreeName << " in "
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
    // create a new tree for holding sorted detector events
    TTree* outputDetectorTree = new TTree(detectorTreeName.c_str(),"");
    if(!outputDetectorTree)
    {
        cerr << "Error: couldn't create detector tree " << detectorTreeName << " when attempting to assign events to macropulses. Exiting... " << endl;
        return 1;
    }
    cout << detectorTreeName << " opened successfully." << endl;
    /**************************************************************************/

    /**************************************************************************/
    // create struct for holding event data and link it to trees
    DetectorEvent detectorEvent;
    
    macropulseTree->SetBranchAddress("macroNo",&detectorEvent.macroNo);
    macropulseTree->SetBranchAddress("macroTime",&detectorEvent.macroTime);
    macropulseTree->SetBranchAddress("targetPos",&detectorEvent.targetPos);

    inputDetectorTree->SetBranchAddress("timetag",&detectorEvent.timetag);
    inputDetectorTree->SetBranchAddress("extTime",&detectorEvent.extTime);
    inputDetectorTree->SetBranchAddress("fineTime",&detectorEvent.fineTime);
    inputDetectorTree->SetBranchAddress("sgQ",&detectorEvent.sgQ);
    inputDetectorTree->SetBranchAddress("lgQ",&detectorEvent.lgQ);
    inputDetectorTree->SetBranchAddress("baseline",&detectorEvent.baseline);
    inputDetectorTree->SetBranchAddress("waveform",&detectorEvent.waveform);

    outputDetectorTree->Branch("macroTime",&detectorEvent.macroTime,"macroTime/D");
    outputDetectorTree->Branch("macroNo",&detectorEvent.macroNo,"macroNo/i");
    outputDetectorTree->Branch("targetPos",&detectorEvent.targetPos,"targetPos/i");

    outputDetectorTree->Branch("completeTime",&detectorEvent.completeTime,"completeTime/d");
    outputDetectorTree->Branch("timetag",&detectorEvent.timetag,"timetag/d");
    outputDetectorTree->Branch("extTime",&detectorEvent.extTime,"extTime/d");
    outputDetectorTree->Branch("fineTime",&detectorEvent.fineTime,"fineTime/d");
    outputDetectorTree->Branch("eventNo",&detectorEvent.eventNo,"eventNo/i");
    outputDetectorTree->Branch("sgQ",&detectorEvent.sgQ,"sgQ/i");
    outputDetectorTree->Branch("lgQ",&detectorEvent.lgQ,"lgQ/i");
    outputDetectorTree->Branch("waveform",&detectorEvent.waveform);
    /**************************************************************************/

    long evtNo = 0;

    long detectorTreeEntries = inputDetectorTree->GetEntries();
    long macropulseEntries = macropulseTree->GetEntries();

    long unsigned int i = 0;
    long unsigned int j = 0;

    inputDetectorTree->GetEntry(j);
    determineDetectorEventTime(detectorEvent);

    double prevMacroTime = 0;

    while(i<macropulseEntries)
    {
        macropulseTree->GetEntry(i);
        detectorEvent.eventNo = 0;

        if(detectorEvent.macroTime < prevMacroTime)
        {
            // mode change
            double prevCompleteTime = 0;

            if(detectorEvent.completeTime > detectorEvent.macroTime+experimentalConfig.facilityConfig.MACRO_LENGTH)
            {
                while(prevCompleteTime < detectorEvent.completeTime)
                {
                    prevCompleteTime = detectorEvent.completeTime;
                    j++;

                    if(j>=detectorTreeEntries)
                    {
                        cout << "Reached end of detector tree. Stopping event assignment to macropulses." << endl;
                        outputDetectorTree->Write();
                        outputFile->Close();
                        inputFile->Close();
                        return 0;
                    }

                    inputDetectorTree->GetEntry(j);
                    determineDetectorEventTime(detectorEvent);
                }
            }
        }

        while(detectorEvent.completeTime < detectorEvent.macroTime)
        {
            j++;

            if(j>=detectorTreeEntries)
            {
                cout << "Reached end of detector tree. Stopping event assignment to macropulses." << endl;
                outputDetectorTree->Write();
                outputFile->Close();
                inputFile->Close();
                return 0;
            }

            inputDetectorTree->GetEntry(j);

            determineDetectorEventTime(detectorEvent);
        }

        while(detectorEvent.completeTime >= detectorEvent.macroTime
           && detectorEvent.completeTime < detectorEvent.macroTime + experimentalConfig.facilityConfig.MACRO_LENGTH)
        {
            outputDetectorTree->Fill();
            detectorEvent.eventNo++;

            j++;

            if(j>=detectorTreeEntries)
            {
                cout << "Reached end of detector tree. Stopping event assignment to macropulses." << endl;
                outputDetectorTree->Write();
                outputFile->Close();
                inputFile->Close();
                return 0;
            }

            inputDetectorTree->GetEntry(j);

            determineDetectorEventTime(detectorEvent);
        }

        // check for complete time wrap-around
        if(detectorEvent.completeTime < detectorEvent.macroTime)
        {
            double prevMacroTime = 0;

            while(detectorEvent.macroTime > prevMacroTime)
            {
                prevMacroTime = detectorEvent.macroTime;

                i++;
                if(i>=macropulseEntries)
                {
                    cout << "Reached end of macropulse tree. Stopping event assignment to macropulses." << endl;
                    outputDetectorTree->Write();
                    outputFile->Close();
                    inputFile->Close();
                    return 0;
                }

                macropulseTree->GetEntry(i);
            }

            continue;
        }

        i++;

        if(i%100==0)
        {
            cout << "processed macropulse " << i << " events.\r";
            fflush(stdout);
        }

        prevMacroTime = detectorEvent.macroTime;
    }

    outputDetectorTree->Write();
    outputFile->Close();
    inputFile->Close();
    return 0;
}
