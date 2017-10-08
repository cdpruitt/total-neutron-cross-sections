#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TDirectoryFile.h"
#include "TROOT.h"

#include "../include/physicalConstants.h"
#include "../include/dataStructures.h"
#include "../include/branches.h"
#include "../include/plots.h"
#include "../include/fillCSHistos.h"
#include "../include/waveform.h"
#include "../include/config.h"
#include "../include/GammaCorrection.h"

using namespace std;

extern Config config;

int fillCSHistos(string inputFileName, string treeName, string outputFileName)
{
    cout << "Filling advanced histograms for tree \"" << treeName << "\"..." << endl;

    // calculate gamma correction for each macropulse
    vector<GammaCorrection> gammaCorrectionList;
    calculateGammaCorrection(inputFileName, treeName, gammaCorrectionList);

    // create TOF histos used for cross section calculation
    vector<TH1D*> TOFHistos;

    for(string targetName : config.targetConfig.TARGET_ORDER)
    {
        string TOFName = targetName + "TOF";
        TOFHistos.push_back(new TH1D(TOFName.c_str(),
                    TOFName.c_str(),
                    config.plotConfig.TOF_BINS,
                    config.plotConfig.TOF_LOWER_BOUND,
                    config.plotConfig.TOF_UPPER_BOUND));
    }

    // create other diagnostic histograms used to examine run data
    TH2D* triangle = new TH2D("triangle","TOF vs. lgQ",config.plotConfig.TOF_RANGE,
            config.plotConfig.TOF_LOWER_BOUND,config.plotConfig.TOF_UPPER_BOUND,4096,0,65536);
    TH1D* timeDiffHisto = new TH1D("time since last event","time since last event",
            config.plotConfig.TOF_RANGE,0,config.plotConfig.TOF_RANGE);
    TH2D* timeDiffVEnergy1 = new TH2D("time difference vs. energy of first",
            "time difference vs. energy of first",config.plotConfig.TOF_RANGE,
            0,config.plotConfig.TOF_RANGE,500,2,700);

    TH2D* time1Vtime2 = new TH2D("time of first vs. time of second",
            "time of first vs. time of second",config.plotConfig.TOF_RANGE,0,
            config.plotConfig.TOF_RANGE,config.plotConfig.TOF_RANGE,0,
            config.plotConfig.TOF_RANGE);

    TH2D* energy1VEnergy2 = new TH2D("energy of first vs. energy of second",
            "energy of first vs. energy of second",
            5*config.plotConfig.NUMBER_ENERGY_BINS, floor(config.plotConfig.ENERGY_LOWER_BOUND), ceil(config.plotConfig.ENERGY_UPPER_BOUND),
            5*config.plotConfig.NUMBER_ENERGY_BINS, floor(config.plotConfig.ENERGY_LOWER_BOUND), ceil(config.plotConfig.ENERGY_UPPER_BOUND));

    TH1D *microNoH = new TH1D("microNoH","microNo",360,0,360);
    TH1D *gammaOffsetH = new TH1D("gammaOffsetH","gammaOffsetH",1000,-3,3);
    TH1D* numberOfGammasH = new TH1D("numberOfGammasH", "number of gammas in each macropulse", 20, 0, 20);
    TH2D* gammaOffsetByGammaNumberH = new TH2D("gammaOffsetByGammaNumber","gammaOffsetByGammaNumber",20,0,20,1000,-5,5);

    TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
    if(!inputFile->IsOpen())
    {
        cerr << "Error: failed to open " << inputFileName << "  to fill histos." << endl;
        return 1;
    }

    TTree* tree = (TTree*)inputFile->Get(treeName.c_str());
    if(!tree)
    {
        cerr << "Error: tried to populate advanced histos, but failed to find " << treeName << " in " << inputFileName << endl;
        return 1;
    }

    DetectorEvent event;

    tree->SetBranchAddress("macroTime",&event.macroTime);
    tree->SetBranchAddress("completeTime",&event.completeTime);
    tree->SetBranchAddress("macroNo",&event.macroNo);
    tree->SetBranchAddress("targetPos",&event.targetPos);
    tree->SetBranchAddress("lgQ",&event.lgQ);

    double prevCompleteTime = 0;
    double prevlgQ = 0;

    double microTime;
    double prevMicroTime = 0;
    int microNo;

    double timeDiff;
    double eventTimeDiff = 0;
    double velocity;
    double rKE;
    double prevRKE = 0;

    unsigned int totalEntries = tree->GetEntries();

    for(long i=0; i<totalEntries; i++)
    {
        tree->GetEntry(i);

        // throw away events during target changer movement
        if(event.targetPos==0)
        {
            continue;
        }

        /*****************************************************************/
        // Calculate event properties

        // find which micropulse the event is in and the time since the start of
        // the micropulse (the TOF)
        timeDiff = event.completeTime-(event.macroTime+gammaCorrectionList[event.macroNo].averageGammaOffset);
        eventTimeDiff = event.completeTime-prevCompleteTime;
        microNo = floor(timeDiff/config.facilityConfig.MICRO_LENGTH);
        microTime = fmod(timeDiff,config.facilityConfig.MICRO_LENGTH);

        // convert micropulse time into neutron velocity based on flight path distance
        velocity = (pow(10.,7.)*config.facilityConfig.FLIGHT_DISTANCE)/microTime; // in meters/sec 

        // convert velocity to relativistic kinetic energy
        rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

        TOFHistos[event.targetPos-1]->Fill(microTime);

        // fill histograms with event data
        triangle->Fill(microTime,event.lgQ);
        timeDiffHisto->Fill(eventTimeDiff);
        timeDiffVEnergy1->Fill(eventTimeDiff,prevRKE);
        time1Vtime2->Fill(prevMicroTime,microTime);
        energy1VEnergy2->Fill(prevRKE,rKE);
        microNoH->Fill(microNo);

        prevlgQ = event.lgQ;
        prevMicroTime = microTime;
        prevCompleteTime = event.completeTime;
        prevRKE = rKE;

        if(i%10000==0)
        {
            cout << "Processed " << i << " " << treeName << " events into CS histos...\r";
            //fflush(stdout);
        }
    }

    cout << endl << "Finished populating \"" << treeName << "\" events into CS histos." << endl;
    cout << "Total events processed = " << totalEntries << endl;

    // write histograms in output file
    TFile* outputFile = new TFile(outputFileName.c_str(),"UPDATE");
    TDirectory* directory = (TDirectory*)outputFile->GetDirectory(treeName.c_str());
    if(!directory)
    {
        directory = outputFile->mkdir(treeName.c_str(),treeName.c_str());
    }

    outputFile->cd(treeName.c_str());

    for(auto histo : TOFHistos)
    {
        histo->Write();
    }

    triangle->Write();
    timeDiffHisto->Write();
    timeDiffVEnergy1->Write();
    time1Vtime2->Write();
    energy1VEnergy2->Write();

    microNoH->Write();

    // record gamma correction data in histograms
    cout << "Adding " << gammaCorrectionList.size() << " gamma corrections to histograms..." << endl;
    for(GammaCorrection gc : gammaCorrectionList)
    {
        gammaOffsetH->Fill(gc.averageGammaOffset);
        numberOfGammasH->Fill(gc.numberOfGammas);
        gammaOffsetByGammaNumberH->Fill(gc.numberOfGammas,gc.averageGammaOffset);
    }

    cout << "Done." << endl << endl;

    gammaOffsetH->Write();
    numberOfGammasH->Write();
    gammaOffsetByGammaNumberH->Write();

    outputFile->Close();

    return 0;
}
