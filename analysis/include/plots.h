#ifndef PLOTS_H
#define PLOTS_H

#include <string>

#include "TH1D.h"
#include "TFile.h"

class Plots
{
    public:
        Plots(std::string name);
        Plots(std::string name, TFile*& inputFile, std::string directory);

        TH1D* getTOFHisto();
        TH1D* getRawTOFHisto();
        TH1D* getEnergyHisto();
        TH1D* getDeadtimeHisto();

    private:
        TH1D* TOFHisto;
        TH1D* rawTOFHisto;
        TH1D* energyHisto;
        TH1D* deadtimeHisto;
};

TH1D* timeBinsToRKEBins(TH1D *inputHisto, std::string name);

#endif
