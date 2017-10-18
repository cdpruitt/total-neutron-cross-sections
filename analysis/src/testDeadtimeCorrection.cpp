#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TDirectoryFile.h"

#include "../include/correctForDeadtime.h"
#include "../include/plots.h"

using namespace std;

int main()
{
    TFile* outputFile = new TFile("testDeadtimeCorrection.root","RECREATE");

    // create model histogram
    TH1D* model = new TH1D("model", "model", 1000, 0, 1000);
    unsigned int numberOfBins = model->GetNbinsX();

    for(int i=450; i<=numberOfBins-450; i++)
    {
        model->SetBinContent(i+1,100);
    }

    model->Write();

    TH1D* measured = new TH1D("measured", "measured", 1000, 0, 1000);

    double numberOfPeriods = 500;
    double deadtimeBins = 0;
    double deadtimeTransitionBins = 0;

    vector<double> fractionDead(numberOfBins);

    for(int i=0; i<numberOfBins; i++)
    {
        for(int j=i-(deadtimeBins+deadtimeTransitionBins); j<i; j++)
        {
            int k=i-(deadtimeBins+deadtimeTransitionBins);

            if((j-k)<deadtimeTransitionBins)
            {
                // deadtime transition region
                if(j<0)
                {
                    fractionDead[i] += (1-fractionDead[i])*(model->GetBinContent(j+numberOfBins+1)/numberOfPeriods)*((j-k)/(double)deadtimeTransitionBins);
                }

                else
                {
                    fractionDead[i] += (1-fractionDead[i])*(model->GetBinContent(j+1)/numberOfPeriods)*((j-k)/(double)deadtimeTransitionBins);
                }
            }

            else
            {
                if(j<0)
                {
                    fractionDead[i] += (1-fractionDead[i])*(model->GetBinContent(j+numberOfBins+1)/numberOfPeriods);
                }

                else
                {
                    fractionDead[i] += (1-fractionDead[i])*(model->GetBinContent(j+1)/numberOfPeriods);
                }
            }
        }

        fractionDead[i] += ((model->GetBinContent(i+1)/numberOfPeriods)/2)*(1-fractionDead[i]);
    }

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
