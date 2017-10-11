#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

#include "../include/config.h"
#include "../include/GammaCorrection.h"
#include "../include/dataStructures.h"
#include "../include/physicalConstants.h"

using namespace std;

extern Config config;

int calculateGammaCorrection(string inputFileName, string treeName, vector<GammaCorrection>& gammaCorrectionList)
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

    const double GAMMA_TIME = pow(10,7)*config.facilityConfig.FLIGHT_DISTANCE/C;
    cout << "For flight distance of " << config.facilityConfig.FLIGHT_DISTANCE
        << " cm, gamma time is " << GAMMA_TIME << "." << endl;

    unsigned int long totalEntries = tree->GetEntries();
    tree->GetEntry(totalEntries-1);
    gammaCorrectionList.resize(macroNo+1, GammaCorrection());

    // Define the range of times considered to be gamma rays
    double timeDiff;
    double microTime;

    for(long i=0; i<totalEntries; i++)
    {
        tree->GetEntry(i);

        timeDiff = completeTime-macroTime;
        microTime = fmod(timeDiff,config.facilityConfig.MICRO_LENGTH);

        // test if gamma
        if(fabs(microTime-GAMMA_TIME)
                <(config.timeOffsetsConfig.GAMMA_WINDOW_SIZE/2))
        {
            gammaCorrectionList[macroNo].averageGammaTime += microTime;
            gammaCorrectionList[macroNo].numberOfGammas++;
        }

        if(i%10000==0)
        {
            cout << "Processed " << i << " events through gamma correction calculation...\r";
            fflush(stdout);
        }
    }

    unsigned int numberOfAverages = 0;
    double overallAverageGammaTime = 0;

    for(int i=0; i<gammaCorrectionList.size(); i++)
    {
        // calculate average gamma offset for each macropulse
        if(gammaCorrectionList[i].numberOfGammas==0)
        {
            continue;
        }

        gammaCorrectionList[i].averageGammaTime /=
             gammaCorrectionList[i].numberOfGammas;

        numberOfAverages++;
        overallAverageGammaTime += gammaCorrectionList[i].averageGammaTime;
    }

    if(numberOfAverages>0)
    {
        overallAverageGammaTime /= (double)numberOfAverages;
    }

    // for all macropulses where no gammas were found to calculate a
    // correction, set the correction to the value of the average correction
    // of all other macropulses
    for(int i=0; i<gammaCorrectionList.size(); i++)
    {
        if(gammaCorrectionList[i].numberOfGammas==0)
        {
            gammaCorrectionList[i].averageGammaTime = overallAverageGammaTime;
        }
    }

    cout << "Finished gamma average calculation."
         << "Gamma average calculated for " << (100*numberOfAverages)/(double)gammaCorrectionList.size()
         << "% of macros." << endl;
    
    cout << "Overall average gamma time was " << overallAverageGammaTime << " ns." << endl;
    cout << "Overall average gamma time - facility gamma time = "
        << overallAverageGammaTime - GAMMA_TIME << " ns." << endl;

    inputFile->Close();

    return 0;
}
