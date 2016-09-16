#ifndef PLOTS_H
#define PLOTS_H

#include <string>

#include "TH1I.h"
#include "TFile.h"

class Plots
{
    public:
        Plots(std::string name);
        Plots(std::string name, TFile*& inputFile);

        TH1I* getTOFHisto();
        TH1I* getEnergyHisto();
        TH1I* getDeadtimeHisto();

    private:
        TH1I* TOFHisto;
        TH1I* energyHisto;
        TH1I* deadtimeHisto;
};

#endif
