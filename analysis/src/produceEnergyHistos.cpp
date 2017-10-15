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

int produceEnergyHistos(string inputFileName, string channelName, string outputFileName)
{
    // open input file
    TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
    if(!inputFile->IsOpen())
    {
        cerr << "Error: failed to open " << inputFileName << "  to fill histos." << endl;
        return 1;
    }

    TDirectory* detectorDirectory = (TDirectory*)inputFile->GetDirectory(channelName.c_str());
    if(!detectorDirectory)
    {
        cerr << "Error: failed to find " << channelName << " directory in " << inputFileName << "." << endl;
        return 1;
    }

    TDirectory* macroTimeDirectory = (TDirectory*)inputFile->GetDirectory("macroTime");
    if(!macroTimeDirectory)
    {
        cerr << "Error: failed to find " << "macroTime" << " directory in " << inputFileName << "." << endl;
        return 1;
    }

    // create output file
    TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");
    TDirectory* outputDirectory = outputFile->mkdir(channelName.c_str(),channelName.c_str());

    // find TOF histograms in input file
    for(int i=0; i<config.targetConfig.TARGET_ORDER.size(); i++)
    {
        string targetName = config.targetConfig.TARGET_ORDER[i];

        string TOFHistoName = targetName + "TOF";
        TH1D* TOF = (TH1D*)detectorDirectory->Get(TOFHistoName.c_str());

        string targetPosHName = "targetPosH";
        TH1I* targetPosH = (TH1I*)macroTimeDirectory->Get(targetPosHName.c_str());

        unsigned int numberOfMacros = targetPosH->GetBinContent(i+2);

        string correctedName = targetName + "Iteration0";
        TH1D* correctedTOF = (TH1D*)TOF->Clone(correctedName.c_str());

        correctedTOF->Write();

        correctForDeadtime(TOF, correctedTOF, numberOfMacros, outputDirectory);

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

