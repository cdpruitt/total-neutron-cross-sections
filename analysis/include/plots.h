#ifndef PLOTS_H
#define PLOTS_H

#include <vector>
#include <string>

#include "TH1I.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "../include/crossSection.h"

class Plots
{
    public:
        Plots(std::string name, TFile*& histoFile, TFile*& waveformFile);

        TH1I* getTOFHisto();
        TH1I* getEnergyHisto();
        TH1I* getTOFHistoCorrected();
        TH1I* getEnergyHistoCorrected();

        TH1I* getWaveformTOFHisto();
        TH1I* getWaveformEnergyHisto();
        TH1I* getWaveformTOFHistoCorrected();
        TH1I* getWaveformEnergyHistoCorrected();

        TH1I* getDeadtimeHisto();

        long getMonitorCounts();
        void setMonitorCounts(int m);

        CrossSection getCrossSection();
        void createCSGraph();

    private:
        TH1I* TOFHisto;
        TH1I* TOFHistoCorrected;
        TH1I* energyHisto;
        TH1I* energyHistoCorrected;
        TGraphErrors* CSGraph;
        TGraphErrors* CSGraphScaledToLit;
        TH1I* deadtimeHisto;

        TH1I* waveformTOFHisto;
        TH1I* waveformTOFHistoCorrected;
        TH1I* waveformEnergyHisto;
        TH1I* waveformEnergyHistoCorrected;
        TGraphErrors* waveformCSGraph;
        TGraphErrors* waveformCSGraphScaledToLit;
        TH1I* waveformDeadtimeHisto;

        // number of counts in monitor during target-in-beam
        // (run-specific)
        long monitorCounts;
        // number of waveform-mode waveforms collected during target-in-beam
        long waveformCounts;

        // calculated cross section
        CrossSection crossSection;
        CrossSection scaledToLit;
};

#endif
