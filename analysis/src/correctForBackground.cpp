#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TDirectory.h"

#include "../include/physicalConstants.h"
#include "../include/correctForBackground.h"
#include "../include/config.h"

using namespace std;

extern Config config;

int correctForBackground(CSPrereqs& csp)
{
    // create corrected histo
    string correctedTOFName = csp.TOFHisto->GetName();
    correctedTOFName = correctedTOFName + "BGCorrected";
    TH1D* correctedTOF = (TH1D*)csp.TOFHisto->Clone(correctedTOFName.c_str());

    int numberOfBins = csp.TOFHisto->GetNbinsX();

    const double GAMMA_TIME = pow(10,7)*config.facility.FLIGHT_DISTANCE/C;
    const double GAMMA_WINDOW_WIDTH = 0.8;

    // "micropulse contamination" correction
    /*
    string gammaName = csp.TOFHisto->GetName();
    gammaName = gammaName + "GammaFit";
    TH1D* gammaFitHisto = (TH1D*)csp.TOFHisto->Clone(gammaName.c_str());

    TF1* gammaFit = new TF1("gammaFit","gaus",(GAMMA_TIME)-GAMMA_WINDOW_WIDTH, (GAMMA_TIME)+GAMMA_WINDOW_WIDTH);
    gammaFitHisto->Fit("gammaFit","Q","",(GAMMA_TIME)-GAMMA_WINDOW_WIDTH, (GAMMA_TIME)+GAMMA_WINDOW_WIDTH);
    double gammaAmplitude = gammaFit->GetParameter(0);

    if(gammaAmplitude <= 0)
    {
        cerr << "Error: failed to fit gamma peak while correcting for background." << endl;
        return 1;
    }

    // fit lobes
    vector<double> offsets = {5, -5};
    vector<TH1D*> lobeHistos;

    for(int j=0; j<offsets.size(); j++)
    {
        double offset = offsets[j];

        string lobeName = csp.target.getName();
        lobeName = lobeName + "Lobe" + to_string(offset);
        TH1D* TOFLobe = (TH1D*)csp.TOFHisto->Clone(lobeName.c_str()); 

        TF1* lobeFit = new TF1("lobeFit","gaus",(GAMMA_TIME+offset)-GAMMA_WINDOW_WIDTH, (GAMMA_TIME+offset)+GAMMA_WINDOW_WIDTH);
        TOFLobe->Fit("lobeFit","Q","",(GAMMA_TIME+offset)-GAMMA_WINDOW_WIDTH, (GAMMA_TIME+offset)+GAMMA_WINDOW_WIDTH);

        double lobeAmplitude = lobeFit->GetParameter(0);
        if(lobeAmplitude <= 0)
        {
            cerr << "Error: failed to fit gamma peak lobe while correcting for background." << endl;
            return 1;
        }

        double ratio = (lobeAmplitude/gammaAmplitude);

        int binOffset = (lobeFit->GetParameter(1)-gammaFit->GetParameter(1))*config.plot.TOF_BINS_PER_NS;

        for(int i=1; i<numberOfBins; i++)
        {
            if(i-binOffset < 0 || i+binOffset>=numberOfBins)
            {
                TOFLobe->SetBinContent(i+binOffset,0);
            }

            TOFLobe->SetBinContent(i+binOffset, csp.TOFHisto->GetBinContent(i)*ratio);
        }

        lobeHistos.push_back(TOFLobe);
    }

    for(auto& lobeHisto : lobeHistos)
    {
        correctedTOF->Add(lobeHisto, -1);
    }
*/
    // "dark current" correction

    // calculate background correction
    double backgroundCounts = 0;
    int backgroundBins = 0;

    for(int i=1; i<=(GAMMA_TIME-10)*config.plot.TOF_BINS_PER_NS; i++)
    {
        backgroundCounts += csp.TOFHisto->GetBinContent(i);
        backgroundBins++;
    }

    backgroundCounts /= backgroundBins;

    for(int i=1; i<=numberOfBins; i++)
    {
        correctedTOF->SetBinContent(
                i, correctedTOF->GetBinContent(i)-backgroundCounts);
    }

    csp.TOFHisto = correctedTOF;

    return 0;
}
