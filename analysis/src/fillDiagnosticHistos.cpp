#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TDirectoryFile.h"
#include "TROOT.h"
#include "TMath.h"

#include "../include/physicalConstants.h"
#include "../include/dataStructures.h"
#include "../include/branches.h"
#include "../include/plots.h"
#include "../include/fillDiagnosticHistos.h"
#include "../include/waveform.h"
#include "../include/config.h"

using namespace std;

extern Config config;

int fillDiagnosticHistos(string inputFileName, string treeName, string outputFileName)
{
    cout << "Filling diagnostic histograms for tree \"" << treeName
        << "\"..." << endl;

    TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
    if(!inputFile->IsOpen())
    {
        cerr << "Error: failed to open " << inputFileName
            << "  to fill histos." << endl;
        return 1;
    }

    TTree* tree = (TTree*)inputFile->Get(treeName.c_str());
    if(!tree)
    {
        cerr << "Error: failed to find " << treeName << " in "
            << inputFileName << endl;
        return 1;
    }

    TFile* outputFile = new TFile(outputFileName.c_str(),"UPDATE");

    TDirectory* directory = (TDirectory*)outputFile->GetDirectory(treeName.c_str());
    if(!directory)
    {
        directory = outputFile->mkdir(treeName.c_str(),treeName.c_str());
    }

    directory->cd();

    // create histos for visualizing basic event data
    TH1D* cycleNumberH = new TH1D("cycleNumberH","cycleNumberH",500,0,500);
    TH1D* macroNoH = new TH1D("macroNoH","macroNo",200000,0,200000);
    TH1D* eventNoH = new TH1D("eventNoH","eventNo",300,0,300);
    TH1D* targetPosH = new TH1D("targetPosH","targetPos",7,0,7);
    TH1D* fineTimeH = new TH1D("fineTimeH","fineTimeH",6200,-2,60);
    TH1D* sgQH = new TH1D("sgQH","sgQ",3500,0,35000);
    TH1D* lgQH = new TH1D("lgQH","lgQ",7000,0,70000);

    TH2D *sgQlgQH = new TH2D("sgQlgQH","short gate Q vs. long gate Q",2048,0,65536,2048,0,65536);
    TH1D *QRatio = new TH1D("QRatio","short gate Q/long gate Q",1000,0,1);

    // create a subdirectory for holding DPP-mode waveform data
    gDirectory->mkdir("waveformsDir","raw DPP waveforms");
    TDirectory* waveformsDir = (TDirectory*)gDirectory->Get("waveformsDir");

    waveformsDir->cd();

    DetectorEvent event;

    tree->SetBranchAddress("cycleNumber",&event.cycleNumber);
    tree->SetBranchAddress("macroNo",&event.macroNo);
    tree->SetBranchAddress("macroTime",&event.macroTime);
    tree->SetBranchAddress("fineTime",&event.fineTime);
    tree->SetBranchAddress("eventNo",&event.eventNo);
    tree->SetBranchAddress("completeTime",&event.completeTime);
    tree->SetBranchAddress("targetPos",&event.targetPos);
    tree->SetBranchAddress("sgQ",&event.sgQ);
    tree->SetBranchAddress("lgQ",&event.lgQ);
    tree->SetBranchAddress("waveform",&event.waveform);

    long totalEntries = tree->GetEntries();

    int prevMacroNo = 0;
    double prevCompleteTime = 0;

    double timeDiff;
    double microTime;

    // loop through the channel-specific tree and populate histos
    for(int i=0; i<totalEntries; i++)
    {
        tree->GetEntry(i);

        cycleNumberH->Fill(event.cycleNumber);
        macroNoH->Fill(event.macroNo);
        targetPosH->Fill(event.targetPos);
        eventNoH->Fill(event.eventNo);
        fineTimeH->Fill(event.fineTime);
        sgQH->Fill(event.sgQ);
        lgQH->Fill(event.lgQ);

        sgQlgQH->Fill(event.sgQ,event.lgQ);
        QRatio->Fill(event.sgQ/(double)event.lgQ);

        if(i%50000==0)
        {
            cout << "Processed " << i << " events through diagnostic histos...\r";
            fflush(stdout);

            stringstream temp;
            temp << "macroNo " << event.macroNo << ", eventNo " << event.eventNo;
            TH1D* waveformH = new TH1D(temp.str().c_str(),temp.str().c_str(),event.waveform->size(),0,event.waveform->size());

            // loop through waveform data and fill histo
            for(int k=0; (size_t)k<event.waveform->size(); k++)
            {
                waveformH->SetBinContent(k,event.waveform->at(k));
            }

            waveformH->Write();
        }
    }

    cout << endl << "Finished populating \"" << treeName << "\" events into diagnostic histos." << endl;
    cout << "Total events processed = " << totalEntries << endl;

    directory->cd();

    cycleNumberH->Write();
    macroNoH->Write();
    eventNoH->Write();
    targetPosH->Write();
    fineTimeH->Write();
    sgQH->Write();
    lgQH->Write();
    sgQlgQH->Write();
    QRatio->Write();

    outputFile->Close();
    inputFile->Close();

    return 0;
}
