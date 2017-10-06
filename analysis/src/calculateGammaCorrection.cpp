#include "TFile.h"
#include "TTree.h"

#include "../include/config.h"
#include "../include/GammaCorrection.h"
#include "../include/dataStructures.h"
#include "../include/physicalConstants.h"

using namespace std;

extern Config config;

int calculateGammaCorrection(string inputFileName, string treeName, vector<GammaCorrection>& gammaCorrectionList)
{
    cout << "Calculating gamma correction for \"" << treeName << "\"..." << endl;

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
    DetectorEvent event;

    tree->SetBranchAddress("macroNo",&event.macroNo);
    tree->SetBranchAddress("macroTime",&event.macroTime);
    tree->SetBranchAddress("completeTime",&event.completeTime);

    unsigned int long totalEntries = tree->GetEntries();
    tree->GetEntry(totalEntries-1);
    gammaCorrectionList.resize(event.macroNo);

    // Define the range of times considered to be gamma rays
    const double GAMMA_TIME = pow(10,7)*config.facilityConfig.FLIGHT_DISTANCE/C;
    cout << "For flight distance of " << config.facilityConfig.FLIGHT_DISTANCE
         << ", gamma time is " << GAMMA_TIME << "." << endl;

    double timeDiff;
    double microTime;

    for(long i=0; i<totalEntries; i++)
    {
        tree->GetEntry(i);

        timeDiff = event.completeTime-event.macroTime;
        microTime = fmod(timeDiff,config.facilityConfig.MICRO_LENGTH);

        // test if gamma
        if(abs(microTime-GAMMA_TIME)<(config.timeOffsetsConfig.GAMMA_WINDOW_SIZE/2))
        {
            gammaCorrectionList[event.macroNo].numberOfGammas++;
            gammaCorrectionList[event.macroNo].averageGammaOffset += microTime;
        }

        if(i%10000==0)
        {
            cout << "Processed " << i << " events through gamma correction calculation...\r";
        }
    }

    unsigned int numberOfOffsets = 0;
    double overallAverageOffset = 0;

    for(GammaCorrection gc : gammaCorrectionList)
    {
        // calculate average gamma offset for each macropulse
        if(gc.numberOfGammas==0)
        {
            gc.averageGammaOffset = 0;
            continue;
        }

        gc.averageGammaOffset = (gc.averageGammaOffset/(double)gc.numberOfGammas)-GAMMA_TIME;

        numberOfOffsets++;
        overallAverageOffset += gc.averageGammaOffset;
    }

    overallAverageOffset /= numberOfOffsets;

    cout << "Finished gamma correction calculation."
         << "Gamma correction calculated for " << (100*numberOfOffsets)/(double)gammaCorrectionList.size()
         << "% of macros." << endl;
    
    cout << "Overall average correction was " << overallAverageOffset << " ns." << endl;

    inputFile->Close();

    return 0;
}
