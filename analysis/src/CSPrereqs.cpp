#include <iostream>
#include <string>
#include <utility>
#include "TFile.h"
#include "TH1D.h"
#include "TDirectory.h"

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

    string TOFHistoName = name + "Corrected";
    TOFHisto = ((TH1D*)gDirectory->Get(TOFHistoName.c_str()));
}

// for histos
void CSPrereqs::getMonitorCounts(TFile* histoFile, string directory, int targetPosition)
{
    TDirectory* dir = (TDirectory*)histoFile->Get(directory.c_str());
    monitorCounts = ((TH1D*)dir->Get("targetPosH"))->GetBinContent(targetPosition+1);
}

void CSPrereqs::getGoodMacroRatio(TFile* histoFile, string directory, string goodMacroHistoName, string macroHistoName)
{
    TDirectory* dir = (TDirectory*)histoFile->Get(directory.c_str());
    TH1D* goodMacrosH = (TH1D*)dir->Get(goodMacroHistoName.c_str());
    unsigned int numberOfBins = goodMacrosH->GetNbinsX();

    unsigned int goodMacroCount = 0;

    for(int j=1; j<=numberOfBins; j++)
    {
        if(goodMacrosH->GetBinContent(j)>0)
        {
            goodMacroCount++;
        }
    }

    TH1D* macroHisto = (TH1D*)dir->Get(macroHistoName.c_str());
    numberOfBins = macroHisto->GetNbinsX();

    unsigned int macroCount = 0;

    for(int j=1; j<numberOfBins; j++)
    {
        if(macroHisto->GetBinContent(j)>0)
        {
            macroCount++;
        }
    }

    goodMacroRatio = goodMacroCount/(double)macroCount;
}

// readData for histos
void CSPrereqs::readEnergyData(TFile* histoFile, string directory, int targetPosition)
{
    // Find deadtime-corrected energy histo for this target
    string histoName = config.targetConfig.TARGET_ORDER[targetPosition];
    getHisto(histoFile, directory, histoName);
}

void CSPrereqs::readMonitorData(TFile* histoFile, string monitorDirectory, string macropulseDirectory, int targetPosition)
{
    // Find monitor histo for this target
    getMonitorCounts(histoFile, monitorDirectory, targetPosition);

    // Scale monitor counts by number of good macropulses for this target
    string goodMacroHistoName = config.targetConfig.TARGET_ORDER[targetPosition] + "GoodMacros";
    string macroHistoName = config.targetConfig.TARGET_ORDER[targetPosition] + "MacroNumber";

    getGoodMacroRatio(histoFile, macropulseDirectory, goodMacroHistoName, macroHistoName);

    monitorCounts *= goodMacroRatio;
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
    augend.energyHisto->Add(addend.energyHisto);
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
    goodMacroRatio = 0;
}

CSPrereqs::CSPrereqs(string targetDataLocation) : target(targetDataLocation)
{
    TOFHisto = new TH1D("","",config.plotConfig.TOF_BINS,config.plotConfig.TOF_LOWER_BOUND,config.plotConfig.TOF_RANGE),target.getName();
    TOFHisto->SetDirectory(0);
    energyHisto = timeBinsToRKEBins(TOFHisto,target.getName());
    energyHisto->SetDirectory(0);

    monitorCounts = 0;
    goodMacroRatio = 0;
}
