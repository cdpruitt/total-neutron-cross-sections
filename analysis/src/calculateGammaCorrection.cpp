#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"

#include <fstream>

#include "../include/config.h"
#include "../include/GammaCorrection.h"
#include "../include/dataStructures.h"
#include "../include/physicalConstants.h"

using namespace std;

extern Config config;

int calculateGammaCorrection(string inputFileName, ofstream& logFile, string treeName, string outputFileName)
{
    // test if output file already exists
    ifstream f(outputFileName);
    if(f.good())
    {
        cout << outputFileName << " already exists; skipping gamma correction calculation." << endl;
        logFile << outputFileName << " already exists; skipping gamma correction calculation." << endl;
        return 2;
    }

    cout << endl << "Start generating gamma correction for each macropulse..." << endl;
    logFile << endl << "*** Gamma Correction ***" << endl;

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
    const double GAMMA_WINDOW_WIDTH = config.time.GAMMA_WINDOW_SIZE/2;

    // find gamma range
    TH1D* uncorrectedTOF = new TH1D(
            "uncorrectedTOF",
            "uncorrectedTOF",
            config.plot.TOF_BINS,
            config.plot.TOF_LOWER_BOUND,
            config.plot.TOF_UPPER_BOUND);

    unsigned int long totalEntries = tree->GetEntries();
    tree->GetEntry(totalEntries-1);

    double timeDiff;
    double microTime;

    for(long i=0; i<totalEntries; i++)
    {
        tree->GetEntry(i);

        timeDiff = completeTime-macroTime;
        microTime = fmod(timeDiff,config.facility.MICRO_LENGTH);

        uncorrectedTOF->Fill(microTime);

        if(i%10000==0)
        {
            cout << "Processed " << i << " events into uncorrectedTOF...\r";
            fflush(stdout);
        }
    }

    vector<double> sumsOfNeighborhood(config.plot.TOF_BINS,0);

    int binNeighborhood = config.plot.TOF_BINS_PER_NS;

    for(int i=1; i<config.plot.TOF_BINS; i++)
    {
        for(int j=i-binNeighborhood; j<i+binNeighborhood; j++)
        {
            sumsOfNeighborhood[i] += uncorrectedTOF->GetBinContent(j);
        }
    }

    int maxSum = 0;
    int maxBin = 0;

    for(int i=0; i<sumsOfNeighborhood.size(); i++)
    {
        if(sumsOfNeighborhood[i] > maxSum)
        {
            maxSum = sumsOfNeighborhood[i];
            maxBin = i;
        }
    }

    double gammaWindowCenter = ((double)maxBin)/config.plot.TOF_BINS_PER_NS;

    tree->GetEntry(totalEntries-1);
    vector<GammaCorrection> gammaCorrectionList(macroNo+1);

    // create advanced histos
    TH1D* gammaHisto = new TH1D("gamma histo", "gamma histo",
            500, gammaWindowCenter-GAMMA_WINDOW_WIDTH,
            gammaWindowCenter+GAMMA_WINDOW_WIDTH);

    TH2D* gammaHisto2D = new TH2D("gamma histo 2D", "gamma histo 2D",
            50, gammaWindowCenter-GAMMA_WINDOW_WIDTH,
            gammaWindowCenter+GAMMA_WINDOW_WIDTH,
            50, gammaWindowCenter-GAMMA_WINDOW_WIDTH,
            gammaWindowCenter+GAMMA_WINDOW_WIDTH);

    TH1D* gammaHistoDiff = new TH1D("gamma histo diff", "gamma histo diff",
            500, -GAMMA_WINDOW_WIDTH,
            GAMMA_WINDOW_WIDTH);

    TH1D *gammaAverageH = new TH1D("gammaAverageH","gammaAverageH",
            120,gammaWindowCenter-GAMMA_WINDOW_WIDTH,
            gammaWindowCenter-GAMMA_WINDOW_WIDTH);

    TH1D* numberOfGammasH = new TH1D("numberOfGammasH",
            "number of gammas in each macropulse", 35, 0, 35);
    TH2D* gammaAverageByGammaNumberH = new TH2D("gammaAverageByGammaNumber",
            "gammaAverageByGammaNumber",40,0,40,60,gammaWindowCenter-3,gammaWindowCenter+3);

    TH1D* timeAutocorrelation;

    TH1D* gammaAverageDiff = new TH1D("gamma average diff", "gamma average diff",
            100, -GAMMA_WINDOW_WIDTH, GAMMA_WINDOW_WIDTH);

    TH2D* gammaAverageDiffByGammaNumber = new TH2D("gamma average diff, 2D",
            "gamma average diff, 2D", 100, -GAMMA_WINDOW_WIDTH,
            GAMMA_WINDOW_WIDTH, 30, 0 ,30);

    TH1D *gammaMicroNoH = new TH1D("gamma microNoH","gammaMicroNoH",360,0,360);

    TH1D *gammaCorrectionH = new TH1D("gammaCorrection", "gammaCorrection",
            gammaCorrectionList.size(), 0, gammaCorrectionList.size());

    // Define the range of times considered to be gamma rays
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
        if(abs(microTime-gammaWindowCenter)<(GAMMA_WINDOW_WIDTH))
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

    logFile << "Gamma peak center determined as: " << gammaWindowCenter << " ns." << endl;
    logFile << "Calculated gamma average: " << overallAverageGammaTime << " ns." << endl;
    logFile << "Used flight distance of " << config.facility.FLIGHT_DISTANCE
        << " cm to generate expected gamma time of " << GAMMA_TIME << " ns." << endl;
    logFile << "Difference between calculated and expected: " << overallAverageGammaTime - GAMMA_TIME << " ns." << endl;
    logFile << (100*numberOfAverages)/(double)gammaCorrectionList.size()
        << "% of macros were used to calculate an average." << endl;

    for(unsigned int i=0; i<gammaCorrectionList.size(); i++)
    {
        GammaCorrection gc = gammaCorrectionList[i];

        gammaAverageH->Fill(gc.averageGammaTime);
        numberOfGammasH->Fill(gc.gammaList.size());
        gammaAverageByGammaNumberH->Fill(gc.gammaList.size(), gc.averageGammaTime);

        gammaCorrectionH->SetBinContent(i+1, gc.correction);
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

    gammaCorrectionH->Write();

    outputFile->Close();
    inputFile->Close();

    logFile << "*** Finished Gamma Correction ***" << endl;

    return 0;
}
