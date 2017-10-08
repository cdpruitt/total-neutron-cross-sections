#include <iostream>
#include <string>
#include <utility>
#include "TFile.h"
#include "TH1I.h"

#include "../include/target.h"
#include "../include/plots.h"
#include "../include/CSPrereqs.h"
#include "../include/config.h"
#include "../include/experiment.h"

using namespace std;

extern Config config;

void CSPrereqs::getHisto(TFile* histoFile, string directory, string name)
{
    // for histos
    string energyHistoName = name + "Energy";

    // for waveforms
    //string energyHistoName = name + "Energy";

    histoFile->cd(directory.c_str());
    energyHisto = ((TH1D*)gDirectory->Get(energyHistoName.c_str()));

    string TOFHistoName = name + "TOF";
    TOFHisto = ((TH1D*)gDirectory->Get(TOFHistoName.c_str()));
}

// for histos
void CSPrereqs::getMonitorCounts(TFile* histoFile, string directory, int targetPosition)
{
    histoFile->cd(directory.c_str());
    monitorCounts = ((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(targetPosition+2);
}

// for waveforms
void CSPrereqs::getMonitorCounts(string monitorFileName, string directory, int targetPosition)
{
    TFile* monitorFile = new TFile(monitorFileName.c_str(),"READ");
    monitorFile->cd(directory.c_str());
    monitorCounts = ((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(targetPosition+2);
    monitorFile->Close();
}

// readData for histos
void CSPrereqs::readEnergyData(TFile* histoFile, string directory, int targetPosition)
{
    // Find deadtime-corrected energy histo for this target
    string histoName = config.targetConfig.TARGET_ORDER[targetPosition];
    getHisto(histoFile, directory, histoName);
}

void CSPrereqs::readMonitorData(TFile* histoFile, string directory, int targetPosition)
{
    // Find monitor histo for this target
    string histoName = config.targetConfig.TARGET_ORDER[targetPosition];
    getMonitorCounts(histoFile, directory, targetPosition);
}

CSPrereqs operator+(CSPrereqs& augend, CSPrereqs& addend)
{
    if(augend.target.getName() != addend.target.getName())
    {
        cerr << "Error: tried to added CSPrereqs of different targets." << endl;
        exit(1);
    }

    augend.energyHisto->Add(addend.energyHisto);
    augend.TOFHisto->Add(addend.TOFHisto);
    augend.monitorCounts += addend.monitorCounts;

    return augend;
}

CSPrereqs::CSPrereqs(Target t)
{
    target = t;
    TOFHisto = new TH1D("","",config.plotConfig.TOF_BINS,config.plotConfig.TOF_LOWER_BOUND,config.plotConfig.TOF_UPPER_BOUND),target.getName();
    TOFHisto->SetDirectory(0);
    energyHisto = timeBinsToRKEBins(TOFHisto,target.getName());
    energyHisto->SetDirectory(0);
    monitorCounts = 0;
}

CSPrereqs::CSPrereqs(string targetDataLocation) : target(targetDataLocation)
{
    TOFHisto = new TH1D("","",config.plotConfig.TOF_BINS,config.plotConfig.TOF_LOWER_BOUND,config.plotConfig.TOF_RANGE),target.getName();
    TOFHisto->SetDirectory(0);
    energyHisto = timeBinsToRKEBins(TOFHisto,target.getName());
    energyHisto->SetDirectory(0);

    monitorCounts = 0;
}
