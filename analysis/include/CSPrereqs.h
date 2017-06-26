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

        void readData(TFile* histoFile, std::string directory, int targetPosition);
        void readData(TFile* histoFile, std::string directory, int targetPosition, std::string monitorFileName);
void getHisto(TFile* histoFile, std::string directory, std::string name);
        void getMonitorCounts(TFile* histoFile, std::string directory, int targetPosition);
        void getMonitorCounts(std::string monitorFileName, std::string directory, int targetPosition);
        void getTargetData(std::string expName, std::string targetName);
        friend CSPrereqs operator+(CSPrereqs& augend, CSPrereqs& addend);

        Target target;     // physical data for this target
        long monitorCounts;// target-specific counts on monitor for a subrun
        TH1D* energyHisto; // target-specific energy histo, corrected for deadtime
        TH1D* TOFHisto; // target-specific energy histo, corrected for deadtime
};



#endif
