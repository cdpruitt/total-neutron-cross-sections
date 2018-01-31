#include <string>
#include <iostream>

#include "TFile.h"
#include "TH1.h"

using namespace std;

const double RATE = 0.01;

int main()
{
    string inputFileName = "/data1/analysis/15/0000/histos.root";
    string inputFile2Name = "/data1/analysis/15/0000/energy.root";

    TFile* inputFile = new TFile(inputFileName.c_str());
    TFile* inputFile2 = new TFile(inputFile2Name.c_str());

    if(!(inputFile->IsOpen()))
    {
        cerr << "Error: failed to open " << inputFileName << "  to fill histos." << endl;
        return 1;
    }

    if(!(inputFile2->IsOpen()))
    {
        cerr << "Error: failed to open " << inputFile2Name << "  to fill histos." << endl;
        return 1;
    }

    string directoryName = "summedDet";
    TDirectory* directory = (TDirectory*)inputFile->Get(directoryName.c_str());
    if(!directory)
    {
        cerr << "Error: failed to open " << directoryName << "  to access histos." << endl;
        return 1;
    }

    TDirectory* directory2 = (TDirectory*)inputFile2->Get(directoryName.c_str());
    if(!directory2)
    {
        cerr << "Error: failed to open " << directoryName << "  to access histos." << endl;
        return 1;
    }

    string goodMacrosHistoName = "blankGoodMacros";
    TH1I* goodMacrosH = (TH1I*)directory->Get(goodMacrosHistoName.c_str());
    if(!goodMacrosH)
    {
        cerr << "Error: failed to open " << goodMacrosHistoName << "  to create rate histo." << endl;
        return 1;
    }

    directory2->cd();

    string TOFHistoName = "blankCorrected";
    TH1D* TOFHisto = (TH1D*)directory2->Get(TOFHistoName.c_str());
    if(!TOFHisto)
    {
        cerr << "Error: failed to open " << TOFHistoName << "  to create rate histo." << endl;
        return 1;
    }

    int numberOfBins = goodMacrosH->GetNbinsX();
    int numberOfMacros = 0;

    for(int j=1; j<=numberOfBins; j++)
    {
        if(goodMacrosH->GetBinContent(j)>0)
        {
            numberOfMacros++;
        }
    }

    double numberOfMicros = numberOfMacros*250;

    int tofBins = TOFHisto->GetNbinsX();

    string outputFileName = "rateHisto.root";
    string rateHistoName = "rateHisto";

    int numberOfRateHistoBins = 1000;

    TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");
    TH1D* rateHisto = new TH1D(rateHistoName.c_str(), rateHistoName.c_str(),
            numberOfRateHistoBins, 0, numberOfRateHistoBins);

    for(int i=0; i<numberOfRateHistoBins; i++)
    {
        rateHisto->SetBinContent(i+1,RATE);
    }

    rateHisto->Write();

    string outputHistoName = "blankRate";
    TH1D* outputHisto = new TH1D(outputHistoName.c_str(), outputHistoName.c_str(),
            tofBins, 0, tofBins);

    for(int i=0; i<tofBins; i++)
    {
        outputHisto->SetBinContent(i+1,TOFHisto->GetBinContent(i+1)/numberOfMicros);
    }

    cout << "Number of micros used for rate generation = " << numberOfMicros << endl;

    outputHisto->Write();
    outputFile->Close();

    inputFile->Close();

    return 0;
}
