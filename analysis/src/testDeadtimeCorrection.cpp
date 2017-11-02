#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TDirectoryFile.h"

using namespace std;

int main()
{
    TFile* outputFile = new TFile("simulation.root","RECREATE");

    // create model histogram
    /*TH1D* model = new TH1D("model", "model", 1000, 0, 1000);
    unsigned int numberOfBins = model->GetNbinsX();

    for(int i=450; i<=numberOfBins-450; i++)
    {
        model->SetBinContent(i+1,100);
    }

    model->Write();
    */

    TH1D* model = (TH1D*)outputFile->Get("allLiveEvents");
    if(!model)
    {
        cerr << "Error: failed to find allLiveEvents histogram in simulation.root." << endl;
        return 1;
    }

    unsigned int numberOfBins = model->GetNbinsX();

    TH1D* measured = new TH1D("measured", "measured", numberOfBins, 0, numberOfBins);

    double numberOfPeriods = 500;
    double deadtimeBins = 0;
    double deadtimeTransitionBins = 0;

    for(int i=0; i<=numberOfBins; i++)
    {
        measured->SetBinContent(i+1,model->GetBinContent(i+1)*(1-fractionDead[i]));
    }

    TH1D* fractionDeadH = new TH1D("fraction dead", "fraction dead", numberOfBins, 0, numberOfBins);
    for(int i=0; i<fractionDead.size(); i++)
    {
        fractionDeadH->SetBinContent(i+1,fractionDead[i]);
    }

    fractionDeadH->Write();

    measured->Write();

    TH1D* corrected = (TH1D*)measured->Clone("Corrected");
    correctForDeadtimeBob(measured, corrected, deadtimeBins, deadtimeTransitionBins, numberOfPeriods);

    corrected->Write();

    outputFile->Close();

    return 0;
}
