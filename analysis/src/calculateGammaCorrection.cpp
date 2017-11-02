#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"

#include "../include/config.h"
#include "../include/GammaCorrection.h"
#include "../include/dataStructures.h"
#include "../include/physicalConstants.h"

using namespace std;

extern Config config;

int calculateGammaCorrection(string inputFileName, string treeName, vector<GammaCorrection>& gammaCorrectionList, string outputFileName)
{
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

    // re-attach to the channel-specific tree for reading out data
    unsigned int macroNo;
    double macroTime;
    double completeTime;
    unsigned int lgQ;

    tree->SetBranchAddress("macroTime",&macroTime);
    tree->SetBranchAddress("completeTime",&completeTime);
    tree->SetBranchAddress("macroNo",&macroNo);
    tree->SetBranchAddress("lgQ",&lgQ);

    // create outputFile
    TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");
    TDirectory* directory = outputFile->GetDirectory(treeName.c_str());

    if(!directory)
    {
        directory = outputFile->mkdir(treeName.c_str(),treeName.c_str());
    }

    directory->cd();


    const double GAMMA_TIME = pow(10,7)*config.facility.FLIGHT_DISTANCE/C;
    cout << "For flight distance of " << config.facility.FLIGHT_DISTANCE
        << " cm, gamma time is " << GAMMA_TIME << "." << endl;
    const double GAMMA_WINDOW_WIDTH = config.timeOffsets.GAMMA_WINDOW_SIZE/2;

    // create advanced histos
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

    TH1D *gammaMicroNoH = new TH1D("gamma microNoH","gammaMicroNoH",360,0,360);

    unsigned int long totalEntries = tree->GetEntries();
    tree->GetEntry(totalEntries-1);
    gammaCorrectionList.resize(macroNo+1, GammaCorrection());

    // Define the range of times considered to be gamma rays
    double timeDiff;
    double microTime;
    int microNo;

    double weight;
    double energy;

    double prevGammaTime = 0;

    for(long i=0; i<totalEntries; i++)
    {
        tree->GetEntry(i);

        timeDiff = completeTime-macroTime;
        microTime = fmod(timeDiff,config.facility.MICRO_LENGTH);

        // test if gamma:
        // if so, use for correction and populate gamma-specific histos
        if(abs(microTime-GAMMA_TIME)<(GAMMA_WINDOW_WIDTH))
        {
            // weight each gamma by the inverse of the FWHM of the gamma peak of
            // just its energy
            if(lgQ<=0)
            {
                lgQ = 1;
            }

            energy = lgQ/(double)200;
            weight = 1/(0.340626 + 1.03717/(energy) + 3.25392/(energy*energy));
            //weight = 1;

            gammaCorrectionList[macroNo].gammaList
                .push_back(GammaEvent(microTime, energy, weight));

            gammaHisto->Fill(microTime);
            gammaHisto2D->Fill(microTime, prevGammaTime);
            gammaHistoDiff->Fill(microTime-prevGammaTime);
            gammaMicroNoH->Fill(microNo);

            prevGammaTime = microTime;
        }

        if(i%10000==0)
        {
            cout << "Processed " << i << " events through gamma correction calculation...\r";
            fflush(stdout);
        }
    }

    unsigned int numberOfAverages = 0;
    double overallAverageGammaTime = 0;

    for(auto& gc : gammaCorrectionList)
    {
        // calculate average gamma offset for each macropulse
        if(gc.gammaList.size()==0)
        {
            continue;
        }

        gc.averageGammaTime = 0;
        double totalWeight = 0;

        for(auto& gammaEvent : gc.gammaList)
        {
            gc.averageGammaTime += gammaEvent.weight*gammaEvent.time;
            totalWeight += gammaEvent.weight;

            gc.numberOfGammas++;
        }

        gc.averageGammaTime /= totalWeight;

        numberOfAverages++;
        overallAverageGammaTime += gc.averageGammaTime;
    }

    if(numberOfAverages>0)
    {
        overallAverageGammaTime /= (double)numberOfAverages;
    }

    // for all macropulses where no gammas were found to calculate a
    // correction, set the correction to the value of the average correction
    // of all other macropulses
    for(auto& gc : gammaCorrectionList)
    {
        if(gc.numberOfGammas==0)
        {
            gc.averageGammaTime = overallAverageGammaTime;
        }

        gc.correction = gc.averageGammaTime-GAMMA_TIME;
    }

    cerr << "Finished gamma average calculation."
        << "Gamma average calculated for " << (100*numberOfAverages)/(double)gammaCorrectionList.size()
        << "% of macros." << endl;

    cerr << "Overall average gamma time was " << overallAverageGammaTime << " ns." << endl;
    cerr << "Overall average gamma time - facility gamma time = "
        << overallAverageGammaTime - GAMMA_TIME << " ns." << endl;

    for(auto& gc : gammaCorrectionList)
    {
        gammaAverageH->Fill(gc.averageGammaTime);
        numberOfGammasH->Fill(gc.gammaList.size());
        gammaAverageByGammaNumberH->Fill(gc.gammaList.size(), gc.averageGammaTime);
    }

    // fill gamma time correction autocorrelation histogram
    vector<double> autocorrelationBins;
    timeAutocorrelation = new TH1D("time autocorrelation",
            "time autocorrelation", (gammaCorrectionList.size()/1000)-1,
            0, ceil(gammaCorrectionList.size()/(double)100));

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
        for(unsigned int i=0; i+delay<gammaCorrectionList.size(); i++)
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

    gammaMicroNoH->Write();

    gammaHisto->Write();
    gammaHisto2D->Write();
    gammaHistoDiff->Write();

    gammaAverageH->Write();
    numberOfGammasH->Write();
    gammaAverageByGammaNumberH->Write();

    timeAutocorrelation->Write();

    gammaAverageDiff->Write();
    gammaAverageDiffByGammaNumber->Write();

    outputFile->Close();
    inputFile->Close();

    return 0;
}
