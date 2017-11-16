#ifndef CS_PREREQS
#define CS_PREREQS

#include <string>
#include "TFile.h"
#include "TH1I.h"
#include "target.h"

/******************************************************************************/
/* Necessary information to calculate a cross section (CS prerequisites) */

class CSPrereqs
{
    public:
        CSPrereqs(Target t);
        CSPrereqs(std::string targetDataLocation);

        void readEnergyData(TFile* histoFile, std::string directory, int targetPosition);
        void readMonitorData(TFile* histoFile, std::string monitorDirectory, int targetPosition);
        void readMacroData(TFile* macroFile, TFile* monitorFile, std::string directory, int targetPosition);

        void getHisto(TFile* histoFile, std::string directory, std::string name);
        void getMonitorCounts(TFile* histoFile, std::string directory, int targetPosition);
        void getMonitorCounts(std::string monitorFileName, std::string directory, int targetPosition);
        void getTargetData(std::string expName, std::string targetName);
        void getGoodMacroRatio(TFile* histoFile, TFile* monitorFile, std::string directory, std::string goodMacroHistoName, std::string macroHistoName);

        friend CSPrereqs operator+(CSPrereqs& augend, CSPrereqs& addend);

        Target target;     // physical data for this target
        double monitorCounts;// target-specific counts on monitor for a subrun
        double goodMacroRatio;

        TH1D* energyHisto; // target-specific energy histo, corrected for deadtime
        TH1D* TOFHisto; // target-specific energy histo, corrected for deadtime
};



#endif
