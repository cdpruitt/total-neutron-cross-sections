#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TDirectory.h"

#include "../include/physicalConstants.h"
#include "../include/correctForBackground.h"
#include "../include/config.h"

using namespace std;

extern Config config;

int correctForBackground(string inputFileName, ofstream& logFile)
{
    cout << "Applying background correction..." << endl;

    // open input file
    TFile* inputFile = new TFile(inputFileName.c_str(),"UPDATE");
    if(!inputFile->IsOpen())
    {
        cerr << "Error: failed to open " << inputFileName << "  to fill histos." << endl;
        return 1;
    }

    for(auto& channelName : config.cs.DETECTOR_NAMES)
    {
        TDirectory* detectorDirectory = inputFile->GetDirectory(channelName.c_str());
        if(!detectorDirectory)
        {
            cerr << "Error: failed to find " << channelName << " directory in " << inputFileName << "." << endl;
            return 1;
        }

        detectorDirectory->cd();

        // find TOF histograms in input file
        for(int i=0; i<config.target.TARGET_ORDER.size(); i++)
        {
            // open uncorrected histo
            string targetName = config.target.TARGET_ORDER[i];

            string TOFHistoName = targetName + "TOF";
            TH1D* TOF = (TH1D*)detectorDirectory->Get(TOFHistoName.c_str());

            if(!TOF)
            {
                cerr << "Error: failed to find " << TOFHistoName << " in "
                    << channelName << " of " << inputFileName << "." << endl;

                inputFile->Close();
                return 1;
            }

            // create corrected histo
            string correctedTOFName = TOFHistoName + "Corrected";
            TH1D* correctedTOF = (TH1D*)TOF->Clone(correctedTOFName.c_str());

            // generate background correction
            int numberOfBins = TOF->GetNbinsX();
            double backgroundCounts = 0;
            int backgroundBins = 0;

            for(int i=1; i<=50*config.plot.TOF_BINS_PER_NS; i++)
            {
                backgroundCounts += TOF->GetBinContent(i);
                backgroundBins++;
            }

            backgroundCounts /= backgroundBins;

            for(int i=1; i<=numberOfBins; i++)
            {
                correctedTOF->SetBinContent(i,
                        correctedTOF->GetBinContent(i)-backgroundCounts);
            }

            correctedTOF->Write();
        }
    }

    inputFile->Close();

    return 0;
}
