#include "TFile.h"
#include "TTree.h"

#include "../include/config.h"
#include "../include/GammaCorrection.h"
#include "../include/dataStructures.h"
#include "../include/physicalConstants.h"

using namespace std;

extern Config config;

int calculateGammaCorrection(string inputFileName, string treeName,
        vector<GammaCorrection>& gammaCorrectionList)
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

    tree->SetBranchAddress("macroTime",&macroTime);
    tree->SetBranchAddress("completeTime",&completeTime);
    tree->SetBranchAddress("macroNo",&macroNo);

    unsigned int long totalEntries = tree->GetEntries();
    tree->GetEntry(totalEntries-1);
    gammaCorrectionList.resize(macroNo);

    // Define the range of times considered to be gamma rays
    const double GAMMA_TIME = pow(10,7)*config.facilityConfig.FLIGHT_DISTANCE/C;
    cout << "For flight distance of " << config.facilityConfig.FLIGHT_DISTANCE
         << ", gamma time is " << GAMMA_TIME << "." << endl;

    double timeDiff;
    double microTime;

    for(long i=0; i<totalEntries; i++)
    {
        tree->GetEntry(i);

        timeDiff = completeTime-macroTime;
        microTime = fmod(timeDiff,config.facilityConfig.MICRO_LENGTH);

        // test if gamma
        if(abs(microTime-GAMMA_TIME)<(config.timeOffsetsConfig.GAMMA_WINDOW_SIZE/2))
        {
            gammaCorrectionList[macroNo].numberOfGammas++;
            gammaCorrectionList[macroNo].averageGammaOffset += microTime;
        }

        if(i%10000==0)
        {
            cout << "Processed " << i << " events through gamma correction calculation...\r";
            fflush(stdout);
        }
    }

    unsigned int numberOfOffsets = 0;
    double overallAverageOffset = 0;

    for(int i=0; i<gammaCorrectionList.size(); i++)
    {
        // calculate average gamma offset for each macropulse
        if(gammaCorrectionList[i].numberOfGammas==0)
        {
            gammaCorrectionList[i].averageGammaOffset = 0;
            continue;
        }

        gammaCorrectionList[i].averageGammaOffset =
            (gammaCorrectionList[i].averageGammaOffset
             /(double)(gammaCorrectionList[i].numberOfGammas))-GAMMA_TIME;

        numberOfOffsets++;
        overallAverageOffset += gammaCorrectionList[i].averageGammaOffset;
    }

    if(numberOfOffsets>0)
    {
        overallAverageOffset /= numberOfOffsets;
    }

    cout << "Finished gamma correction calculation."
         << "Gamma correction calculated for " << (100*numberOfOffsets)/(double)gammaCorrectionList.size()
         << "% of macros." << endl;
    
    cout << "Overall average correction was " << overallAverageOffset << " ns." << endl;

    inputFile->Close();

    return 0;
}
