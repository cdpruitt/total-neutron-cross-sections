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
#include "TRandom3.h"

#include "../include/physicalConstants.h"
#include "../include/dataStructures.h"
#include "../include/branches.h"
#include "../include/plots.h"
#include "../include/fillBasicHistos.h"
#include "../include/waveform.h"
#include "../include/config.h"

using namespace std;

extern ProcessedEvent procEvent;

extern Config config;

int fillBasicHistos(string inputFileName, string treeName, string outputFileName)
{
    cout << "Filling basic histograms for tree \"" << treeName << "\"..." << endl;

    TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
    if(!inputFile->IsOpen())
    {
        cerr << "Error: failed to open " << inputFileName << "  to fill histos." << endl;
        return 1;
    }

    TTree* tree = (TTree*)inputFile->Get(treeName.c_str());
    if(!tree)
    {
        cerr << "Error: tried to populate basic histos, but failed to find " << treeName << " in " << inputFileName << endl;
        return 1;
    }

    TFile* outputFile = new TFile(outputFileName.c_str(),"UPDATE");

    outputFile->mkdir(treeName.c_str(),treeName.c_str());
    outputFile->GetDirectory(treeName.c_str())->cd();

    // instantiate DPP-mode histograms
    TH1I* macroNoH = new TH1I("macroNoH","macroNo",200000,0,200000);
    macroNoH->GetXaxis()->SetTitle("macropulse number of each event");

    TH1I* macroNoHBlank = new TH1I("macroNoHBlank","macroNoBlank",200000,0,200000);
    macroNoHBlank->GetXaxis()->SetTitle("macropulse number of each event, blank target");

    TH1I* macroNoHTarget1 = new TH1I("macroNoHTarget1","macroNoTarget1",200000,0,200000);
    macroNoHTarget1->GetXaxis()->SetTitle("macropulse number of each event, target 1");

    TH1I* macroNoHTarget2 = new TH1I("macroNoHTarget2","macroNoTarget2",200000,0,200000);
    macroNoHTarget2->GetXaxis()->SetTitle("macropulse number of each event, target 2");

    TH1I* macroNoHTarget3 = new TH1I("macroNoHTarget3","macroNoTarget3",200000,0,200000);
    macroNoHTarget3->GetXaxis()->SetTitle("macropulse number of each event, target 3");

    TH1I* macroNoHTarget4 = new TH1I("macroNoHTarget4","macroNoTarget4",200000,0,200000);
    macroNoHTarget4->GetXaxis()->SetTitle("macropulse number of each event, target 4");

    TH1I* macroNoHTarget5 = new TH1I("macroNoHTarget5","macroNoTarget5",200000,0,200000);
    macroNoHTarget5->GetXaxis()->SetTitle("macropulse number of each event, target 5");

    TH1I* eventNoH = new TH1I("eventNoH","eventNo",300,0,300);
    eventNoH->GetXaxis()->SetTitle("event number of each event");

    TH1I* macroTimeH = new TH1I("macroTimeH","macroTime",6000,0,6000000000);
    macroTimeH->GetXaxis()->SetTitle("macropulse time zero for each event");

    TH1I* targetPosH = new TH1I("targetPosH","targetPos",7,0,7);
    targetPosH->GetXaxis()->SetTitle("target position of each event");

    TH1I* sgQH = new TH1I("sgQH","sgQ",3500,0,35000);
    sgQH->GetXaxis()->SetTitle("short gate integrated charge for each event");

    TH1I* lgQH = new TH1I("lgQH","lgQ",7000,0,70000);
    lgQH->GetXaxis()->SetTitle("long gate integrated charge for each event");

    TH1I* macroNoDiffH = new TH1I("macroNoDiffH","macroNoDiff",100,0,100);
    macroNoDiffH->GetXaxis()->SetTitle("difference between macropulse numbers of consecutive events");

    TH1I* completeTimeH = new TH1I("completeTimeH","completeTime",pow(2,20),0,pow(2,32));
    completeTimeH->GetXaxis()->SetTitle("complete time of event");

    TH1I* fineTimeH = new TH1I("fineTimeH","fineTimeH",6200,-2,60);
    fineTimeH->GetXaxis()->SetTitle("fine time of event");

    TH1I* diffCompleteTimeH = new TH1I("diffCompleteTimeH","diffCompleteTime",pow(2,16),0,pow(2,26));
    diffCompleteTimeH->GetXaxis()->SetTitle("difference between complete time of consecutive events");

    TH1I* diffMacroCompleteTimesH = new TH1I("diffMacroCompleteTimesH","diffMacroCompleteTime",pow(2,20),0,pow(2,20));
    diffMacroCompleteTimesH->GetXaxis()->SetTitle("difference between complete time of event and its macrotime");

    vector<TH1D*> TOFHistos;

    for(unsigned int i=0; i<config.targetConfig.TARGET_ORDER.size(); i++)
    {
        string TOFName = config.targetConfig.TARGET_ORDER[i] + "TOFBasic";
        TOFHistos.push_back(new TH1D(TOFName.c_str(),TOFName.c_str(),config.plotConfig.TOF_BINS,config.plotConfig.TOF_LOWER_BOUND,config.plotConfig.TOF_UPPER_BOUND));
    }

    // create a subdirectory for holding DPP-mode waveform data
    gDirectory->mkdir("waveformsDir","raw DPP waveforms");
    TDirectory* waveformsDir = (TDirectory*)gDirectory->Get("waveformsDir");
    waveformsDir->cd();

    ProcessedEvent procEvent;
    setBranchesHistos(tree, procEvent);

    long totalEntries = tree->GetEntries();

    cout << "Populating " << treeName << " histograms..." << endl;

    int prevMacroNo = 0;
    double prevCompleteTime = 0;

    double timeDiff;
    double microTime;

    // loop through the channel-specific tree and populate histos
    for(long j=0; j<totalEntries; j++)
    {
        tree->GetEntry(j);

        timeDiff = procEvent.completeTime-procEvent.macroTime;
        microTime = fmod(timeDiff,config.facilityConfig.MICRO_LENGTH);

        if(j%50000==0)
        {
            cout << "Processed " << j << " events through basic histos...\r";
            fflush(stdout);

            stringstream temp;
            temp << "macroNo " << procEvent.macroNo << ", eventNo " << procEvent.eventNo;
            TH1I* waveformH = new TH1I(temp.str().c_str(),temp.str().c_str(),procEvent.waveform->size(),0,procEvent.waveform->size());

            // loop through waveform data and fill histo
            for(int k=0; (size_t)k<procEvent.waveform->size(); k++)
            {
                waveformH->SetBinContent(k,procEvent.waveform->at(k));
            }

            waveformH->Write();
        }

        switch(procEvent.targetPos)
        {
            case 1:
                TOFHistos[0]->Fill(microTime);
                macroNoHBlank->Fill(procEvent.macroNo);
                break;
            case 2:
                TOFHistos[1]->Fill(microTime);
                macroNoHTarget1->Fill(procEvent.macroNo);
                break;
            case 3:
                TOFHistos[2]->Fill(microTime);
                macroNoHTarget2->Fill(procEvent.macroNo);
                break;
            case 4:
                TOFHistos[3]->Fill(microTime);
                macroNoHTarget3->Fill(procEvent.macroNo);
                break;
            case 5:
                TOFHistos[4]->Fill(microTime);
                macroNoHTarget4->Fill(procEvent.macroNo);
                break;
            case 6:
                TOFHistos[5]->Fill(microTime);
                macroNoHTarget5->Fill(procEvent.macroNo);
                break;
            default:
                break;
        }

        macroNoH->Fill(procEvent.macroNo);
        targetPosH->Fill(procEvent.targetPos);
        macroTimeH->Fill(procEvent.macroTime);

        macroNoDiffH->Fill(procEvent.macroNo-prevMacroNo);

        completeTimeH->Fill(procEvent.completeTime);

        fineTimeH->Fill(procEvent.fineTime);
        diffCompleteTimeH->Fill(procEvent.completeTime-prevCompleteTime);

        diffMacroCompleteTimesH->Fill(procEvent.completeTime-procEvent.macroTime);

        eventNoH->Fill(procEvent.eventNo);
        sgQH->Fill(procEvent.sgQ);
        lgQH->Fill(procEvent.lgQ);

        prevMacroNo = procEvent.macroNo;
        prevCompleteTime = procEvent.completeTime;

    }

    outputFile->Write();

    inputFile->Close();
    outputFile->Close();

    return 0;
}
