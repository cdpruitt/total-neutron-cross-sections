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

int produceEnergyHistos(string inputFileName, ofstream& log, string channelName, string outputFileName)
{
    // open input file
    TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
    if(!inputFile->IsOpen())
    {
        cerr << "Error: failed to open " << inputFileName << "  to fill histos." << endl;
        return 1;
    }

    TDirectory* detectorDirectory = inputFile->GetDirectory(channelName.c_str());
    if(!detectorDirectory)
    {
        cerr << "Error: failed to find " << channelName << " directory in " << inputFileName << "." << endl;
        return 1;
    }

    // create output file
    TFile* outputFile = new TFile(outputFileName.c_str(),"UPDATE");
    TDirectory* outputDirectory = outputFile->GetDirectory(channelName.c_str());
    if(!outputDirectory)
    {
        outputDirectory = outputFile->mkdir(channelName.c_str(),channelName.c_str());
    }

    outputDirectory->cd();

    // find TOF histograms in input file
    for(int i=0; i<config.target.TARGET_ORDER.size(); i++)
    {
        string targetName = config.target.TARGET_ORDER[i];

        string TOFHistoName = targetName + "TOF";
        TH1D* TOF = (TH1D*)detectorDirectory->Get(TOFHistoName.c_str());

        string goodMacrosHistoName = targetName + "GoodMacros";
        TH1I* goodMacrosH = (TH1I*)detectorDirectory->Get(goodMacrosHistoName.c_str());
        unsigned int numberOfBins = goodMacrosH->GetNbinsX();
        unsigned int numberOfMacros = 0;

        for(int j=1; j<=numberOfBins; j++)
        {
            if(goodMacrosH->GetBinContent(j)>0)
            {
                numberOfMacros++;
            }
        }

        double numberOfMicros = numberOfMacros*(config.facility.MICROS_PER_MACRO);

        string correctedName = targetName + "Uncorrected";
        TH1D* correctedTOF = (TH1D*)TOF->Clone(correctedName.c_str());

        correctedTOF->Write();

        outputDirectory->cd();
        correctForDeadtime(TOF, correctedTOF, numberOfMicros);

        correctedName = targetName + "Corrected";
        correctedTOF = (TH1D*)correctedTOF->Clone(correctedName.c_str());

        correctedTOF->Write();

        TH1D* correctedEnergy = convertTOFtoEnergy(correctedTOF, targetName + "Energy");
        correctedEnergy->Write();
    }
    
    inputFile->Close();
    outputFile->Close();

    return 0;
}
