#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TDirectoryFile.h"

#include "../include/config.h"
#include "../include/correctForDeadtime.h"
#include "../include/plots.h"

using namespace std;

extern Config config;

int produceEnergyHistos(CSPrereqs& csp)
{
    cout << endl << "Mapping " << csp.name << " TOF histogram to energy domain..." << endl;

    // open input file
    TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
    if(!inputFile->IsOpen())
    {
        cerr << "Error: failed to open " << inputFileName << "  to fill histos." << endl;
        return 1;
    }

    // create output file
    TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");

    for(auto& channelName : config.cs.DETECTOR_NAMES)
    {
        TDirectory* detectorDirectory = inputFile->GetDirectory(channelName.c_str());
        if(!detectorDirectory)
        {
            cerr << "Error: failed to find " << channelName << " directory in " << inputFileName << "." << endl;
            return 1;
        }

        TDirectory* outputDirectory = outputFile->mkdir(channelName.c_str(),channelName.c_str());

        // find TOF histograms in input file
        for(int i=0; i<config.target.TARGET_ORDER.size(); i++)
        {
            string targetName = config.target.TARGET_ORDER[i];
            string TOFHistoName = targetName + "TOFCorrected";

            detectorDirectory->cd();

            TH1D* TOF = (TH1D*)detectorDirectory->Get(TOFHistoName.c_str());
            if(!TOF)
            {
                cerr << "Error: failed to open " << TOFHistoName << " in " << inputFileName << endl;
                return 1;
            }

            outputDirectory->cd();

            TH1D* correctedEnergy = convertTOFtoEnergy(TOF, targetName + "Energy");
            TOF->Write();
            correctedEnergy->Write();
        }
    }

    outputFile->Close();
    inputFile->Close();

    return 0;
}
