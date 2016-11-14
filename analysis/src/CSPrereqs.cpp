#include <iostream>
#include <string>
#include <utility>
#include "TFile.h"
#include "TH1I.h"

#include "../include/target.h"
#include "../include/CSPrereqs.h"
#include "../include/analysisConstants.h"

using namespace std;

void CSPrereqs::getHisto(TFile* histoFile, string directory, string name)
{
    histoFile->cd(directory.c_str());
    energyHisto = ((TH1I*)gDirectory->Get(name.c_str()));
}

void CSPrereqs::getMonitorCounts(TFile* histoFile, string directory, int targetPosition)
{
    histoFile->cd(directory.c_str());
    monitorCounts = ((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(targetPosition+2);
}

void CSPrereqs::readData(TFile* histoFile, string directory, int targetPosition)
{
    // Find deadtime-corrected energy histo for this target
    string histoName = positionNames[targetPosition] + "CorrectedEnergy";
    getHisto(histoFile, directory, histoName);

    // Find number of events in the monitor for each target to use in scaling
    // cross-sections
    getMonitorCounts(histoFile, "monitor", targetPosition);
}

CSPrereqs operator+(CSPrereqs& augend, CSPrereqs& addend)
{
    if(augend.target.getName() != addend.target.getName())
    {
        cerr << "Error: tried to added CSPrereqs of different targets." << endl;
        exit(1);
    }

    augend.energyHisto->Add(addend.energyHisto);
    augend.monitorCounts += addend.monitorCounts;

    return augend;
}

CSPrereqs::CSPrereqs(Target t)
{
    target = t;
    energyHisto = 0;
    monitorCounts = 0;
}

CSPrereqs::CSPrereqs(string targetDataLocation) : target(targetDataLocation)
{
    energyHisto = 0;
    monitorCounts = 0;
}
