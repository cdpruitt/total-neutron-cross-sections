#include "TFile.h"
#include "TTree.h"

using namespace std;

int produceEnergyHistos(string inputFileName, string channelName, string energyFileName)
{
    // open input file
    TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
    if(!inputFile->IsOpen())
    {
        cerr << "Error: failed to open " << inputFileName << "  to fill histos." << endl;
        return 1;
    }

    TDirectory* directory = (TDirectory*)inputFile->GetDirectory(channelName.c_str());
    if(!directory)
    {
        cerr << "Error: failed to find " << channelName << " directory in " << inputFileName << "." << endl;
        return 1;
    }
    directory->cd();

    // create output file
    TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");
    outputFile->mkdir(channelName.c_str(),channelName.c_str());
    outputFile->cd(treeName.c_str());

    // find TOF histograms in input file
    for(string targetName : config.targetConfig.TARGET_ORDER)
    {
        string TOFhistoName = targetName + "TOF";
        string macroNoHistoName = targetName + "macroNo";
        TH1D* TOF = (TH1D*)gDirectory->Get(histoName.c_str());
        TH1I* macroNoH = (TH1I*)gDirectory->Get(macroNoHName.c_str());
        TH1D* deadtimeCorrectedTOF = generateDeadtimeCorrection(TOF, macroNoH);
        TH1D* correctedEnergy = convertTOFtoEnergy(TOF, targetName + "Energy");
    correctedEnergy->Write();



    return 0;
}

