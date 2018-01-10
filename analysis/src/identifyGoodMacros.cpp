#include "../include/identifyGoodMacros.h"
#include "../include/config.h"

#include "TFile.h"
#include "TTree.h"

#include <iostream>

using namespace std;

extern Config config;

int identifyGoodMacros(string inputFileName, vector<MacropulseEvent>& macropulseList, ofstream& logFile)
{
    ifstream f(inputFileName);
    if(!f.good())
    {
        cout << inputFileName << " does not exist; cannot identify good macros." << endl;
        logFile << inputFileName << " does not exist; cannot identify good macros. " << endl;
        return 2;
    }

    cout << "Starting identification of good macropulses..." << endl;

    TFile* inputFile = new TFile(inputFileName.c_str(), "READ");
    if(!inputFile->IsOpen())
    {
        cerr << "Error: failed to open " << inputFileName << "  to fill histos." << endl;
        return 1;
    }

    TTree* eventTree =
        (TTree*)inputFile->Get(config.analysis.GAMMA_CORRECTION_TREE_NAME.c_str());

    TTree* monitorTree =
        (TTree*)inputFile->Get(config.analysis.MONITOR_TREE_NAME.c_str());

    if(!eventTree || !monitorTree)
    {
        cerr << "Error: couldn't find either event tree or monitor tree when trying to identify good macropulses." << endl;
        inputFile->Close();
        return 1;
    }

    // connect event tree to event data buffer
    DetectorEvent event;
    vector<int>* waveformPointer = 0;

    eventTree->SetBranchAddress("cycleNumber",&event.cycleNumber);
    eventTree->SetBranchAddress("macroNo",&event.macroNo);
    eventTree->SetBranchAddress("macroTime",&event.macroTime);
    eventTree->SetBranchAddress("fineTime",&event.fineTime);
    eventTree->SetBranchAddress("eventNo",&event.eventNo);
    eventTree->SetBranchAddress("completeTime",&event.completeTime);
    eventTree->SetBranchAddress("targetPos",&event.targetPos);
    eventTree->SetBranchAddress("sgQ",&event.sgQ);
    eventTree->SetBranchAddress("lgQ",&event.lgQ);
    eventTree->SetBranchAddress("waveform",&waveformPointer);

    DetectorEvent monitorEvent;
    monitorTree->SetBranchAddress("macroNo",&monitorEvent.macroNo);

    // tally the number of events in each macropulse
    eventTree->GetEntry(0);
    macropulseList.push_back(MacropulseEvent(event.macroNo,0,event.targetPos));

    unsigned int totalEntries = eventTree->GetEntries();
    for(long i=1; i<totalEntries; i++)
    {
        eventTree->GetEntry(i);

        if(macropulseList.back().macroNo<event.macroNo)
        {
            macropulseList.push_back(MacropulseEvent(event.macroNo, 0, event.targetPos));
        }

        macropulseList.back().numberOfEventsInMacro++;

        if(i%10000==0)
        {
            cout << "Found " << macropulseList.size() << " macropulses during good macropulse identification...\r";
            fflush(stdout);
        }
    }

    unsigned int currentEntry = 0;

    unsigned int monitorTotalEntries = monitorTree->GetEntries();
    for(int i=0; i<monitorTotalEntries; i++)
    {
        monitorTree->GetEntry(i);

        while(macropulseList[currentEntry].macroNo < monitorEvent.macroNo)
        {
            currentEntry++;
        }

        macropulseList[currentEntry].numberOfMonitorsInMacro++;

        if(i%10000==0)
        {
            cout << "Found " << macropulseList.size() << " macropulses during good macropulse identification...\r";
            fflush(stdout);
        }
    }

    cout << endl;

    // calculate the average number of events in each macropulse, by target
    vector<double> averageEventsPerMacropulseByTarget(config.target.TARGET_ORDER.size(),0);
    vector<unsigned int> numberOfMacropulsesByTarget(config.target.TARGET_ORDER.size(),0);

    for(auto& macropulse : macropulseList)
    {
        averageEventsPerMacropulseByTarget[macropulse.targetPos]
            += macropulse.numberOfEventsInMacro;
        numberOfMacropulsesByTarget[macropulse.targetPos]++;
    }

    for(int i=0; i<averageEventsPerMacropulseByTarget.size(); i++)
    {
        averageEventsPerMacropulseByTarget[i] /= numberOfMacropulsesByTarget[i];
        logFile << "Target \"" << config.target.TARGET_ORDER[i]
            << "\" had, on average, " << averageEventsPerMacropulseByTarget[i]
            << " events per macropulse." << endl;
    }

    cout << "Finished calculating average events per macropulse."
        << endl;

    /*TH1D* averageRatePerTarget = new TH1D("average rate per target", "average rate per target",
      config.target.TARGET_ORDER.size(), 0, config.target.TARGET_ORDER.size());

    for(int i=0; i<averageEventsPerMacropulseByTarget.size(); i++)
    {
        averageRatePerTarget->SetBinContent(i+1, averageEventsPerMacropulseByTarget[i]);
    }

    averageRatePerTarget->Write();
    */

    // using the average just calculated, identify macros with too few events
    // per macro (indicating readout problem)
    for(auto& macropulse : macropulseList)
    {
        if((macropulse.numberOfEventsInMacro
              >(0.8)*averageEventsPerMacropulseByTarget[macropulse.targetPos])
         && (macropulse.numberOfMonitorsInMacro > 0)
         )
        /*if((macropulse.numberOfEventsInMacro > 0) &
                (macropulse.numberOfMonitorsInMacro > 0))
                */
        {
            // good macro
            macropulse.isGoodMacro = true;
        }

        else
        {
            macropulse.isGoodMacro = false;
        }
    }

    inputFile->Close();

    return 0;
}
