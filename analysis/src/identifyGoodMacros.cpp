#include "../include/identifyGoodMacros.h"
#include "../include/config.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#include <iostream>

using namespace std;

extern Config config;

int identifyGoodMacros(string macropulseFileName, vector<MacropulseEvent>& macropulseList, ofstream& logFile)
{
    if(macropulseList.size()==0)
    {
        cerr << "Error: cannot identify good macros from empty macro list." << endl;
        return 1;
    }

    // check to see if output file already exists; if so, exit
    ifstream f(macropulseFileName);

    if(f.good())
    {
        cout << macropulseFileName << " already exists; skipping event assignment to macropulses." << endl;
        logFile << macropulseFileName << " already exists; skipping event assignment to macropulses." << endl;
        return 2;
    }

    f.close();

    // calculate the average number of events in each macropulse, by target
    vector<double> averageEventsPerMacropulseByTarget(config.target.TARGET_ORDER.size(),0);
    vector<int> numberOfMacropulsesByTarget(config.target.TARGET_ORDER.size(),0);

    for(auto& macropulse : macropulseList)
    {
        averageEventsPerMacropulseByTarget[macropulse.targetPos]
            += macropulse.numberOfEventsInMacro;
        numberOfMacropulsesByTarget[macropulse.targetPos]++;
    }

    for(int i=0; i<averageEventsPerMacropulseByTarget.size(); i++)
    {
        averageEventsPerMacropulseByTarget[i] /= (double)(numberOfMacropulsesByTarget[i]);
        logFile << "Target \"" << config.target.TARGET_ORDER[i]
            << "\" had, on average, " << averageEventsPerMacropulseByTarget[i]
            << " events per macropulse." << endl;
    }

    cout << "Finished calculating average events per macropulse."
        << endl;

    // using the average just calculated, identify macros with too few events
    // per macro (indicating readout problem)
    for(auto& macropulse : macropulseList)
    {
        if((macropulse.numberOfEventsInMacro
              >(0.5)*averageEventsPerMacropulseByTarget[macropulse.targetPos])
         && (macropulse.numberOfMonitorsInMacro > 0)
         )
        {
            // good macro
            macropulse.isGoodMacro = true;
        }

        else
        {
            macropulse.isGoodMacro = false;
        }
    }

    cout << "Finished identifying good macropulses." << endl;

    TFile* outputFile = new TFile(macropulseFileName.c_str(),"CREATE");

    TH1D* cycleNumber = new TH1D("cycleNumber","cycleNumber", 1000, 0, 1000);
    TH1D* macroNumber = new TH1D("macroNo","macroNo", 500000, 0, 500000);
    TH1D* macroTime = new TH1D("macroTime","macroTime", 5000000, 0, 5000000000);

    vector<TH1D*> macroNumberByTargets;
    vector<TH1D*> eventsPerMacroByTargets;
    vector<TH1D*> monitorsPerMacroByTargets;

    for(auto& targetName : config.target.TARGET_ORDER)
    {
        string macroNumberByTargetsName = targetName + "macroNo";
        macroNumberByTargets.push_back(new TH1D(macroNumberByTargetsName.c_str(), macroNumberByTargetsName.c_str(), 500000, 0, 500000));

        string eventsPerMacroName = targetName + "eventsPerMacro";
        eventsPerMacroByTargets.push_back(new TH1D(eventsPerMacroName.c_str(), eventsPerMacroName.c_str(), 300, 0, 300));

        string monitorsPerMacroName = targetName + "monitorsPerMacro";
        monitorsPerMacroByTargets.push_back(new TH1D(monitorsPerMacroName.c_str(), monitorsPerMacroName.c_str(), 100, 0, 100));
    }

    for(auto& macropulse : macropulseList)
    {
        cycleNumber->Fill(macropulse.cycleNumber);
        macroNumber->Fill(macropulse.macroNo);
        macroTime->Fill(macropulse.macroTime);

        macroNumberByTargets[macropulse.targetPos]->Fill(macropulse.macroNo);
        eventsPerMacroByTargets[macropulse.targetPos]->Fill(macropulse.numberOfEventsInMacro);
        monitorsPerMacroByTargets[macropulse.targetPos]->Fill(macropulse.numberOfMonitorsInMacro);
    }

    TTree* macropulseTree = new TTree("macropulses","");

    MacropulseEvent me;

    macropulseTree->Branch("cycleNumber",&me.cycleNumber,"cycleNumber/I");
    macropulseTree->Branch("macroNo",&me.macroNo,"macroNo/I");
    macropulseTree->Branch("macroTime",&me.macroTime,"macroTime/d");
    macropulseTree->Branch("targetPos",&me.targetPos,"targetPos/I");
    macropulseTree->Branch("numberOfEventsInMacro",&me.numberOfEventsInMacro,"numberOfEventsInMacro/I");
    macropulseTree->Branch("numberOfMonitorsInMacro",&me.numberOfMonitorsInMacro,"numberOfMonitorsInMacro/I");
    macropulseTree->Branch("isGoodMacro",&me.isGoodMacro,"isGoodMacro/O");

    for(auto& macropulse : macropulseList)
    {
        me.cycleNumber = macropulse.cycleNumber;
        me.macroNo = macropulse.macroNo;
        me.macroTime = macropulse.macroTime;
        me.targetPos = macropulse.targetPos;
        me.numberOfEventsInMacro = macropulse.numberOfEventsInMacro;
        me.numberOfMonitorsInMacro = macropulse.numberOfMonitorsInMacro;
        me.isGoodMacro = macropulse.isGoodMacro;

        macropulseTree->Fill();
    }

    macropulseTree->Write();

    cycleNumber->Write();
    macroNumber->Write();
    macroTime->Write();

    for(int i=0; i<macroNumberByTargets.size(); i++)
    {
        macroNumberByTargets[i]->Write();
        eventsPerMacroByTargets[i]->Write();
        monitorsPerMacroByTargets[i]->Write();
    }

    outputFile->Close();

    return 0;
}
