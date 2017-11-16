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

    histoFile->cd(directory.c_str());
    energyHisto = ((TH1D*)gDirectory->Get(energyHistoName.c_str()));

    if(!energyHisto)
    {
        cerr << "Error: Failed to open histogram \"" << energyHistoName
            << "\" in file " << histoFile
            << " (in getHisto)." << endl;
        return;
    }

    string TOFHistoName = name + "TOFCorrected";
    TOFHisto = ((TH1D*)gDirectory->Get(TOFHistoName.c_str()));

    if(!TOFHisto)
    {
        cerr << "Error: Failed to open histogram \"" << TOFHistoName
            << "\" in file " << histoFile
            << " (in getHisto)." << endl;
        return;
    }
}

// for histos
void CSPrereqs::getMonitorCounts(TFile* histoFile, string directory, int targetPosition)
{
    TDirectory* dir = (TDirectory*)histoFile->Get(directory.c_str());
    if(!dir)
    {
        cerr << "Error: Failed to open directory " << directory << " in file " << histoFile->GetName()
            << " (in getMonitorCounts)." << endl;
        return;
    }

    TH1D* targetPosHisto = (TH1D*)dir->Get("targetPosH");
    if(!targetPosHisto)
    {
        cerr << "Error: Failed to open histogram \"targetPosH\" in file " << histoFile
            << " (in getMonitorCounts)." << endl;
        return;
    }

    monitorCounts = targetPosHisto->GetBinContent(targetPosition+1);

    if(monitorCounts<=0)
    {
        cerr << "Error: for target position " << config.target.TARGET_ORDER[targetPosition] << ", monitor counts read was <=0 (in getMonitorCounts)." << endl;
        return;
    }
}

void CSPrereqs::getGoodMacroRatio(TFile* histoFile, TFile* monitorFile, string directory, string goodMacroHistoName, string macroHistoName)
{
    TDirectory* dir = (TDirectory*)histoFile->Get(directory.c_str());
    if(!dir)
    {
        cerr << "Error: Failed to open directory " << directory << " in file " << histoFile->GetName()
            << " (in getMonitorCounts)." << endl;
        return;
    }

    TDirectory* monitorDir = (TDirectory*)monitorFile->Get(directory.c_str());
    if(!monitorDir)
    {
        cerr << "Error: Failed to open directory " << directory << " in file " << monitorFile->GetName()
            << " (in getMonitorCounts)." << endl;
        return;
    }

    TH1D* goodMacrosH = (TH1D*)dir->Get(goodMacroHistoName.c_str());
    if(!goodMacrosH)
    {
        cerr << "Error: Failed to open histogram \"" << goodMacroHistoName
            << "\" in file " << histoFile->GetName() << " (in getGoodMacroRatio)." << endl;
        return;
    }

    unsigned int numberOfBins = goodMacrosH->GetNbinsX();

    unsigned int goodMacroCount = 0;

    for(int j=1; j<=numberOfBins; j++)
    {
        if(goodMacrosH->GetBinContent(j)>0)
        {
            goodMacroCount++;
        }
    }

    TH1D* macroHisto = (TH1D*)monitorDir->Get(macroHistoName.c_str());
    if(!macroHisto)
    {
        cerr << "Error: Failed to open histogram \"" << macroHistoName
            << "\" in file " << histoFile->GetName() << " (in getGoodMacroRatio)." << endl;
        return;
    }

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

    if(goodMacroRatio<=0 || goodMacroRatio>=1)
    {
        cerr << "Error: goodMacroRatio was " << goodMacroRatio << ", outside bounds [0,1]." << endl;
    }
}

// readData for histos
void CSPrereqs::readEnergyData(TFile* histoFile, string directory, int targetPosition)
{
    // Find deadtime-corrected energy histo for this target
    string histoName = config.target.TARGET_ORDER[targetPosition];
    getHisto(histoFile, directory, histoName);
}

void CSPrereqs::readMonitorData(TFile* histoFile, string monitorDirectory, int targetPosition)
{
    // Find monitor histo for this target
    getMonitorCounts(histoFile, monitorDirectory, targetPosition);
}

void CSPrereqs::readMacroData(TFile* macroFile, TFile* monitorFile, string detectorName, int targetPosition)
{
    string goodMacroHistoName = config.target.TARGET_ORDER[targetPosition] + "GoodMacros";
    string macroHistoName = config.target.TARGET_ORDER[targetPosition] + "MacroNumber";

    getGoodMacroRatio(macroFile, monitorFile, detectorName, goodMacroHistoName, macroHistoName);
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
    TOFHisto = new TH1D("","",config.plot.TOF_BINS,config.plot.TOF_LOWER_BOUND,config.plot.TOF_UPPER_BOUND),target.getName();
    TOFHisto->SetDirectory(0);
    energyHisto = timeBinsToRKEBins(TOFHisto,target.getName());
    energyHisto->SetDirectory(0);

    monitorCounts = 0;
    goodMacroRatio = 0;
}

CSPrereqs::CSPrereqs(string targetDataLocation) : target(targetDataLocation)
{
    TOFHisto = new TH1D("","",config.plot.TOF_BINS,config.plot.TOF_LOWER_BOUND,config.plot.TOF_RANGE),target.getName();
    TOFHisto->SetDirectory(0);
    energyHisto = timeBinsToRKEBins(TOFHisto,target.getName());
    energyHisto->SetDirectory(0);

    monitorCounts = 0;
    goodMacroRatio = 0;
}
