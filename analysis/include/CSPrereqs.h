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

        void readEnergyData(TFile* histoFile, std::string directory, int targetPosition);
        void readMonitorData(TFile* histoFile, std::string monitorDirectory, int targetPosition);
        void readMacroData(TFile* macroFile, std::string directory, int targetPosition);
        void readEventData(TFile* macroFile, std::string directory, int targetPosition);

        void getHisto(TFile* histoFile, std::string directory, std::string name);
        void getMonitorCounts(TFile* histoFile, std::string directory, int targetPosition);
        void getMonitorCounts(std::string monitorFileName, std::string directory, int targetPosition);
        void getTargetData(std::string expName, std::string targetName);
        void getMacroNumber(TFile* histoFile, std::string directory, std::string goodMacroHistoName, std::string macroHistoName);
        void getAverageRate(TFile* histoFile, std::string averageRateDataName, unsigned int targetNumber);

        friend CSPrereqs operator+(CSPrereqs& augend, CSPrereqs& addend);

        Target target;     // physical data for this target
        double monitorCounts;// target-specific counts on monitor for a subrun
        double goodMacroNumber;
        double totalMacroNumber;
        double totalEventNumber;

        TH1D* energyHisto; // target-specific energy histo, corrected for deadtime
        TH1D* TOFHisto; // target-specific energy histo, corrected for deadtime
};

void extractGraphData(
        TGraphErrors* graph,
        std::vector<double>* xValues,
        std::vector<double>* xError,
        std::vector<double>* yValues,
        std::vector<double>* yError);

void extractGraphData(
        TGraphErrors* graph,
        DataSet& dataSet);

double calculateRMS(std::vector<double> graph1Data, std::vector<double> graph2Data);

DataSet scale(DataSet setToScale, DataSet expReference, DataSet litReference);

int readTargetData(std::vector<CSPrereqs>& allCSPrereqs, std::string expName);

int readSubRun(std::vector<CSPrereqs>& allCSPrereqs, std::string expName, int runNumber, int subRun, std::string detectorName, std::string dataLocation);

#endif
