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

    const double Q_LOW_THRESHOLD = 0;
    const double Q_HIGH_THRESHOLD = 65550;
    const double NUMBER_OF_EVENTS = 50000000;

    // calculate gamma correction for each macropulse
    vector<GammaCorrection> gammaCorrectionList;
    calculateGammaCorrection(inputFileName, treeName, gammaCorrectionList);

    // create outputFile
    TFile* outputFile = new TFile(outputFileName.c_str(),"UPDATE");
    TDirectory* directory = (TDirectory*)outputFile->GetDirectory(treeName.c_str());
    if(!directory)
    {
        directory = outputFile->mkdir(treeName.c_str(),treeName.c_str());
    }
    directory->cd();

    // create histos
    const double GAMMA_TIME = pow(10,7)*config.facilityConfig.FLIGHT_DISTANCE/C;

    TH1D* gammaHisto = new TH1D("gamma histo", "gamma histo",
            500, GAMMA_TIME-config.timeOffsetsConfig.GAMMA_WINDOW_SIZE/2,
            GAMMA_TIME+config.timeOffsetsConfig.GAMMA_WINDOW_SIZE/2);

    TH2D* gammaHisto2D = new TH2D("gamma histo 2D", "gamma histo 2D",
            50, GAMMA_TIME-config.timeOffsetsConfig.GAMMA_WINDOW_SIZE/2,
            GAMMA_TIME+config.timeOffsetsConfig.GAMMA_WINDOW_SIZE/2,
            50, GAMMA_TIME-config.timeOffsetsConfig.GAMMA_WINDOW_SIZE/2,
            GAMMA_TIME+config.timeOffsetsConfig.GAMMA_WINDOW_SIZE/2);

    TH1D* gammaHistoDiff = new TH1D("gamma histo diff", "gamma histo diff",
            500, -config.timeOffsetsConfig.GAMMA_WINDOW_SIZE/2,
            config.timeOffsetsConfig.GAMMA_WINDOW_SIZE/2);

    TH1D *gammaAverageH = new TH1D("gammaAverageH","gammaAverageH",
            120,GAMMA_TIME-config.timeOffsetsConfig.GAMMA_WINDOW_SIZE/2,
            GAMMA_TIME-config.timeOffsetsConfig.GAMMA_WINDOW_SIZE/2);

    TH1D* numberOfGammasH = new TH1D("numberOfGammasH",
            "number of gammas in each macropulse", 35, 0, 35);
    TH2D* gammaAverageByGammaNumberH = new TH2D("gammaAverageByGammaNumber",
            "gammaAverageByGammaNumber",40,0,40,60,GAMMA_TIME-3,GAMMA_TIME+3);

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

    for(auto& gc : gammaCorrectionList)
    {
        if(gc.numberOfGammas==0)
        {
            continue;
        }

        gammaAverageH->Fill(gc.averageGammaTime);
        numberOfGammasH->Fill(gc.numberOfGammas);
        gammaAverageByGammaNumberH->Fill(gc.numberOfGammas, gc.averageGammaTime);
    }

    // fill gamma time correction autocorrelation histogram
    vector<double> autocorrelationBins;
    /*unsigned int i = 0;
    while(pow(1.1,i)<gammaCorrectionList.size())
    {
        autocorrelationBins.push_back(pow(1.1,i));
        i++;
    }

    TH2D* timeAutocorrelation = new TH2D("time autocorrelation",
            "time autocorrelation", autocorrelationBins.size()-1,
            &autocorrelationBins[0],
            config.timeOffsetsConfig.GAMMA_WINDOW_SIZE*40,
            -config.timeOffsetsConfig.GAMMA_WINDOW_SIZE,
            config.timeOffsetsConfig.GAMMA_WINDOW_SIZE);
            */

    unsigned int i = 0;
    TH2D* timeAutocorrelation = new TH2D("time autocorrelation",
            "time autocorrelation", (gammaCorrectionList.size()/1000)-1,
            0, ceil(gammaCorrectionList.size()/(double)100),
            200, -1, 1);

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
                fflush(stdout);
            }
        }

        correlation /= (gammaCorrectionList.size()-1)*gammaAverageVariance;

        timeAutocorrelation->Fill(delay,correlation);
    }

    gammaAverageH->Write();
    numberOfGammasH->Write();
    gammaAverageByGammaNumberH->Write();

    timeAutocorrelation->Write();

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

    double prevGammaTime = 0;

    double prevAverageTime = 0;



    unsigned int totalEntries = tree->GetEntries();

    for(long i=0; i<totalEntries && i<NUMBER_OF_EVENTS; i++)
    {
        tree->GetEntry(i);

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
        if(abs(microTime-GAMMA_TIME)<(config.timeOffsetsConfig.GAMMA_WINDOW_SIZE/2))
        {
            gammaHisto->Fill(microTime);
            gammaHisto2D->Fill(microTime, prevGammaTime);
            gammaHistoDiff->Fill(microTime-prevGammaTime);

            prevGammaTime = microTime;
        }

        TOFHistos[event.targetPos-1]->Fill(microTime);

        // fill histograms with event data
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

        if(i%10000==0)
        {
            cout << "Processed " << i << " " << treeName << " events into CS histos...\r";
            //fflush(stdout);
        }
    }

    cout << endl << "Finished populating \"" << treeName << "\" events into CS histos." << endl;
    cout << "Total events processed = " << totalEntries << endl;

    outputFile->cd(treeName.c_str());

    gammaHisto->Write();
    gammaHisto2D->Write();
    gammaHistoDiff->Write();

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

    outputFile->Close();

    return 0;
}
