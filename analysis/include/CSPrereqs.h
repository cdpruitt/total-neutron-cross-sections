#ifndef CS_PREREQS
#define CS_PREREQS

#include <string>
#include "TFile.h"
#include "TH1I.h"
#include "target.h"
#include "dataSet.h"

/******************************************************************************/
/* Necessary information to calculate a cross section (CS prerequisites) */

class CSPrereqs
{
    public:
        CSPrereqs() {};
        CSPrereqs(Target t);
        CSPrereqs(std::string targetDataLocation);

        int readTOFHisto(TFile* histoFile, std::string directory, std::string targetName);
        int readMonitorCounts(TFile* histoFile, std::string directory, std::string targetName);
        int readMacroData(TFile* macroFile, std::string directory, std::string targetName);
        int readEventData(TFile* macroFile, std::string directory, std::string targetName);
        int readUncorrectedTOFHisto(TFile* macroFile, std::string directory, std::string targetName);

        void readTargetData(std::string expName, std::string targetName);
        void getAverageRate(TFile* histoFile, std::string averageRateDataName, int targetNumber);

        friend CSPrereqs operator+(CSPrereqs& augend, CSPrereqs& addend);

        Target target;     // physical data for this target
        double monitorCounts = 0;// target-specific counts on monitor for a subrun
        double goodMacroNumber = 0;
        double totalMacroNumber = 0;
        double totalEventNumber = 0;

        TH1D* energyHisto = 0; // target-specific energy histo, corrected for deadtime
        TH1D* TOFHisto = 0; // target-specific energy histo, corrected for deadtime
        TH1D* uncorrectedTOFHisto = 0; // target-specific energy histo, uncorrected for deadtime

};

void extractGraphData(
        TGraphAsymmErrors* graph,
        std::vector<double>* xValues,
        std::vector<double>* xError,
        std::vector<double>* yValues,
        std::vector<double>* yError);

void extractGraphData(
        TGraphAsymmErrors* graph,
        DataSet& dataSet);

double calculateRMS(std::vector<double> graph1Data, std::vector<double> graph2Data);

DataSet scale(DataSet setToScale, DataSet expReference, DataSet litReference);

int readTargetData(std::vector<CSPrereqs>& allCSPrereqs, std::string expName);

int readSubRun(CSPrereqs& subRunData, std::string expName, int runNumber, int subRun, std::string detectorName, std::string dataLocation);

#endif
