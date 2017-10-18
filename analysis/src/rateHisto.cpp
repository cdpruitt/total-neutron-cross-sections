#include <string>
#include <iostream>

#include "TFile.h"
#include "TH1.h"

using namespace std;

const double RATE = 0.01;
const double numberOfBins = 1000;

int main(int argc, char** argv)
{
    string outputFileName = argv[1];
    string rateHistoName = argv[2];

    TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");
    TH1D* rateHisto = new TH1D(rateHistoName.c_str(), rateHistoName.c_str(),
            numberOfBins, 0, numberOfBins);

    for(unsigned int i; i<numberOfBins; i++)
    {
        rateHisto->SetBinContent(i+1,RATE);
    }

    rateHisto->Write();

    outputFile->Close();

    return 0;
}
