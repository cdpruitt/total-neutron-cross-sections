#include <iostream>
#include <string>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TDirectoryFile.h"
#include "TRandom3.h"

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

const double Q_LOW_THRESHOLD = 0;
const double Q_HIGH_THRESHOLD = 65550;
const double NUMBER_OF_EVENTS = 50000000;

int fillCSHistos(string inputFileName, string treeName, string outputFileName)
{
    cout << "Filling histograms for tree \"" << treeName << "\"..." << endl;

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

    // connect input tree to event data buffer
    DetectorEvent event;
    vector<int>* waveformPointer = 0;

    tree->SetBranchAddress("cycleNumber",&event.cycleNumber);
    tree->SetBranchAddress("macroNo",&event.macroNo);
    tree->SetBranchAddress("macroTime",&event.macroTime);
    tree->SetBranchAddress("fineTime",&event.fineTime);
    tree->SetBranchAddress("eventNo",&event.eventNo);
    tree->SetBranchAddress("completeTime",&event.completeTime);
    tree->SetBranchAddress("targetPos",&event.targetPos);
    tree->SetBranchAddress("sgQ",&event.sgQ);
    tree->SetBranchAddress("lgQ",&event.lgQ);
    tree->SetBranchAddress("waveform",&waveformPointer);

    // create outputFile
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

    // create histos
    const double GAMMA_TIME = pow(10,7)*config.facilityConfig.FLIGHT_DISTANCE/C;
    const double GAMMA_WINDOW_WIDTH = config.timeOffsetsConfig.GAMMA_WINDOW_SIZE/2;
    TH1D* gammaHisto = new TH1D("gamma histo", "gamma histo",
            500, GAMMA_TIME-GAMMA_WINDOW_WIDTH,
            GAMMA_TIME+GAMMA_WINDOW_WIDTH);

    TH2D* gammaHisto2D = new TH2D("gamma histo 2D", "gamma histo 2D",
            50, GAMMA_TIME-GAMMA_WINDOW_WIDTH,
            GAMMA_TIME+GAMMA_WINDOW_WIDTH,
            50, GAMMA_TIME-GAMMA_WINDOW_WIDTH,
            GAMMA_TIME+GAMMA_WINDOW_WIDTH);

    TH1D* gammaHistoDiff = new TH1D("gamma histo diff", "gamma histo diff",
            500, -GAMMA_WINDOW_WIDTH,
            GAMMA_WINDOW_WIDTH);

    TH1D *gammaAverageH = new TH1D("gammaAverageH","gammaAverageH",
            120,GAMMA_TIME-GAMMA_WINDOW_WIDTH,
            GAMMA_TIME-GAMMA_WINDOW_WIDTH);

    TH1D* numberOfGammasH = new TH1D("numberOfGammasH",
            "number of gammas in each macropulse", 35, 0, 35);
    TH2D* gammaAverageByGammaNumberH = new TH2D("gammaAverageByGammaNumber",
            "gammaAverageByGammaNumber",40,0,40,60,GAMMA_TIME-3,GAMMA_TIME+3);

    TH1D* timeAutocorrelation;

    TH1D* gammaAverageDiff = new TH1D("gamma average diff", "gamma average diff",
            100, -GAMMA_WINDOW_WIDTH, GAMMA_WINDOW_WIDTH);

    TH2D* gammaAverageDiffByGammaNumber = new TH2D("gamma average diff, 2D",
            "gamma average diff, 2D", 100, -GAMMA_WINDOW_WIDTH,
            GAMMA_WINDOW_WIDTH, 30, 0 ,30);

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
            config.plotConfig.TOF_LOWER_BOUND,config.plotConfig.TOF_UPPER_BOUND/5,4096,0,65536);
    TH2D* triangleEnergy = new TH2D("triangleEnergy","Energy vs. lgQ",10*config.plotConfig.NUMBER_ENERGY_BINS,
            config.plotConfig.ENERGY_LOWER_BOUND,config.plotConfig.ENERGY_UPPER_BOUND,4096,0,65536);

    TH1D* timeDiffHisto = new TH1D("time since last event","time since last event",
            config.plotConfig.TOF_RANGE,0,config.plotConfig.TOF_RANGE);
    TH2D* timeDiffVEnergy1 = new TH2D("time difference vs. energy of first",
            "time difference vs. energy of first",config.plotConfig.TOF_RANGE,
            0,config.plotConfig.TOF_RANGE,10*config.plotConfig.NUMBER_ENERGY_BINS,2,700);

    TH2D* time1Vtime2 = new TH2D("time of first vs. time of second",
            "time of first vs. time of second",config.plotConfig.TOF_RANGE,0,
            config.plotConfig.TOF_RANGE,config.plotConfig.TOF_RANGE,0,
            config.plotConfig.TOF_RANGE);

    TH2D* energy1VEnergy2 = new TH2D("energy of first vs. energy of second",
            "energy of first vs. energy of second",
            10*config.plotConfig.NUMBER_ENERGY_BINS, floor(config.plotConfig.ENERGY_LOWER_BOUND), ceil(config.plotConfig.ENERGY_UPPER_BOUND),
            10*config.plotConfig.NUMBER_ENERGY_BINS, floor(config.plotConfig.ENERGY_LOWER_BOUND), ceil(config.plotConfig.ENERGY_UPPER_BOUND));

    TH1D *microNoH = new TH1D("microNoH","microNo",360,0,360);

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

    double prevGammaTime = 0;
    double prevAverageTime = 0;

    // if this channel is a detector
    // channel, we should fill additional detector histograms
    bool fillDetectorHistos = false;

    for(auto& name : config.csConfig.DETECTOR_NAMES)
    {
        if(treeName==name)
        {
            fillDetectorHistos = true;
            break;
        }
    }

    vector<GammaCorrection> gammaCorrectionList;
    if(fillDetectorHistos)
    {
        // calculate gamma correction for each macropulse
        calculateGammaCorrection(inputFileName, treeName, gammaCorrectionList);

        for(auto& gc : gammaCorrectionList)
        {
            gammaAverageH->Fill(gc.averageGammaTime);
            numberOfGammasH->Fill(gc.gammaList.size());
            gammaAverageByGammaNumberH->Fill(gc.gammaList.size(), gc.averageGammaTime);
        }

        // fill gamma time correction autocorrelation histogram
        vector<double> autocorrelationBins;
        /*unsigned int i = 0;
          while(pow(1.1,i)<gammaCorrectionList.size())
          {
          autocorrelationBins.push_back(pow(1.1,i));
          i++;
          }

          TH1D* timeAutocorrelation = new T12D("time autocorrelation",
          "time autocorrelation", autocorrelationBins.size()-1,
          &autocorrelationBins[0]);
          */

        unsigned int i = 0;
        timeAutocorrelation = new TH1D("time autocorrelation",
                "time autocorrelation", (gammaCorrectionList.size()/1000)-1,
                0, ceil(gammaCorrectionList.size()/(double)100));

        // calculate overall average gamma time
        double overallAverageGammaTime = 0;

        for(auto& gc : gammaCorrectionList)
        {
            overallAverageGammaTime += gc.averageGammaTime;
        }

        overallAverageGammaTime /= gammaCorrectionList.size();

        // calculate average gamma time variance
        double gammaAverageVariance = 0;

        for(auto& gc : gammaCorrectionList)
        {
            gammaAverageVariance += pow((gc.averageGammaTime-overallAverageGammaTime),2);
        }

        gammaAverageVariance /= gammaCorrectionList.size();

        // calculate time autocorrelation
        double correlation = 0;

        for(int delay = 100; delay<gammaCorrectionList.size()/50; delay += 100)
        {
            for(i=0; i+delay<gammaCorrectionList.size(); i++)
            {
                correlation +=
                    (gammaCorrectionList[i].averageGammaTime-overallAverageGammaTime)
                    *(gammaCorrectionList[i+delay].averageGammaTime-overallAverageGammaTime);

                if(i%1000==0)
                {
                    cout << "Calculated autocorrelation (delay = " << delay << ") through " << i << " gamma averages...\r";
                }
            }

            correlation /= (gammaCorrectionList.size()-1)*gammaAverageVariance;

            timeAutocorrelation->Fill(delay,correlation);
        }

        TRandom3* rng = new TRandom3();

        // calculate variance of gamma average
        for(GammaCorrection gc : gammaCorrectionList)
        {
            double gammaAverage1 = 0;
            double gammaAverage2 = 0;

            vector<GammaEvent> selectedGammas;

            unsigned int randomGammaNumber = 0;

            while(selectedGammas.size()<gc.gammaList.size())
            {
                randomGammaNumber = floor(rng->Uniform(0, gc.gammaList.size()));
                selectedGammas.push_back(gc.gammaList[randomGammaNumber]);
                gc.gammaList.erase(gc.gammaList.begin()+randomGammaNumber);
            }

            for(auto& ge : selectedGammas)
            {
                gammaAverage1 += ge.time;
            }

            gammaAverage1 /= selectedGammas.size();

            for(auto& ge : gc.gammaList)
            {
                gammaAverage2 += ge.time;
            }

            gammaAverage2 /= gc.gammaList.size();

            gammaAverageDiff->Fill(gammaAverage2-gammaAverage1);
            gammaAverageDiffByGammaNumber->Fill(gammaAverage2-gammaAverage1, gc.gammaList.size());
        }

        cout << endl;
    }

    unsigned int totalEntries = tree->GetEntries();

    for(long i=0; i<totalEntries && i<NUMBER_OF_EVENTS; i++)
    {
        tree->GetEntry(i);

        event.waveform = *waveformPointer;

        cycleNumberH->Fill(event.cycleNumber);
        macroNoH->Fill(event.macroNo);
        targetPosH->Fill(event.targetPos);
        eventNoH->Fill(event.eventNo);
        fineTimeH->Fill(event.fineTime);
        sgQH->Fill(event.sgQ);
        lgQH->Fill(event.lgQ);

        sgQlgQH->Fill(event.sgQ,event.lgQ);
        QRatio->Fill(event.sgQ/(double)event.lgQ);

        if(i%10000==0)
        {
            cout << "Processed " << i << " " << treeName << " events into CS histos...\r";

            waveformsDir->cd();
            stringstream temp;
            temp << "macroNo " << event.macroNo << ", eventNo " << event.eventNo;
            TH1D* waveformH = new TH1D(temp.str().c_str(),temp.str().c_str(),event.waveform.size(),0,event.waveform.size());

            // loop through waveform data and fill histo
            for(int k=0; (size_t)k<event.waveform.size(); k++)
            {
                waveformH->SetBinContent(k,event.waveform[k]);
            }

            waveformH->Write();
        }

        if(!fillDetectorHistos)
        {
            continue;
        }

        // throw away events during target changer movement
        if(event.targetPos==0)
        {
            continue;
        }

        // gates:
        if(event.lgQ<Q_LOW_THRESHOLD || event.lgQ>Q_HIGH_THRESHOLD)
        {
            continue;
        }

        /*****************************************************************/
        // Calculate event properties

        // find which micropulse the event is in and the time since the start of
        // the micropulse (the TOF)
        timeDiff = event.completeTime-event.macroTime;

        // correct times using average gamma time
        timeDiff += (GAMMA_TIME-gammaCorrectionList[event.macroNo].averageGammaTime);

        eventTimeDiff = event.completeTime-prevCompleteTime;
        microNo = floor(timeDiff/config.facilityConfig.MICRO_LENGTH);
        microTime = fmod(timeDiff,config.facilityConfig.MICRO_LENGTH);

        // convert micropulse time into neutron velocity based on flight path distance
        velocity = (pow(10.,7.)*config.facilityConfig.FLIGHT_DISTANCE)/microTime; // in meters/sec 

        // convert velocity to relativistic kinetic energy
        rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

        // test if gamma
        if(abs(microTime-GAMMA_TIME)<(GAMMA_WINDOW_WIDTH))
        {
            gammaHisto->Fill(microTime);
            gammaHisto2D->Fill(microTime, prevGammaTime);
            gammaHistoDiff->Fill(microTime-prevGammaTime);

            prevGammaTime = microTime;
        }

        TOFHistos[event.targetPos-1]->Fill(microTime);

        // fill detector histograms with event data
        triangle->Fill(microTime,event.lgQ);
        triangleEnergy->Fill(rKE,event.lgQ);
        timeDiffHisto->Fill(eventTimeDiff);
        timeDiffVEnergy1->Fill(eventTimeDiff,prevRKE);
        time1Vtime2->Fill(prevMicroTime,microTime);
        energy1VEnergy2->Fill(prevRKE,rKE);
        microNoH->Fill(microNo);

        prevlgQ = event.lgQ;
        prevMicroTime = microTime;
        prevCompleteTime = event.completeTime;
        prevRKE = rKE;
    }

    cout << endl << "Finished populating \"" << treeName << "\" events into CS histos." << endl;
    cout << "Total events processed = " << totalEntries << endl;

    directory->cd();

    if(fillDetectorHistos)
    {
        for(auto histo : TOFHistos)
        {
            histo->Write();
        }

        triangle->Write();
        triangleEnergy->Write();
        timeDiffHisto->Write();
        timeDiffVEnergy1->Write();
        time1Vtime2->Write();
        energy1VEnergy2->Write();
        microNoH->Write();

        gammaHisto->Write();
        gammaHisto2D->Write();
        gammaHistoDiff->Write();

        gammaAverageH->Write();
        numberOfGammasH->Write();
        gammaAverageByGammaNumberH->Write();

        timeAutocorrelation->Write();

        gammaAverageDiff->Write();
        gammaAverageDiffByGammaNumber->Write();
    }

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
