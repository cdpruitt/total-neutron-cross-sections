#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TLatex.h"

#include "../include/target.h"
#include "../include/dataSet.h"
#include "../include/dataPoint.h"
#include "../include/CSPrereqs.h"
#include "../include/CSUtilities.h"
#include "../include/crossSection.h"
#include "../include/experiment.h"
#include "../include/plots.h"

using namespace std;

Config config;

int main(int, char* argv[])
{
    string dataLocation = argv[1];

    string expName = argv[2]; // experiment directory where runs to-be-sorted
                              // are listed

    int runNumber = atoi(argv[3]);

    int lowSubrun = atoi(argv[4]);
    int highSubrun = atoi(argv[5]);

    string detectorName = argv[6]; // detector name to be used for calculating cross sections

    vector<CSPrereqs> allCSPrereqs;
    vector<vector<double>> runningMonitorCounts;
    vector<vector<double>> runningGoodMacroCounts;
 
    // read in run config file
    config = Config(expName, runNumber);

    readTargetData(allCSPrereqs, expName);

    // Loop through all subruns of this run
    for(int subRun=lowSubrun; subRun<=highSubrun; subRun++)
    {
        readSubRun(allCSPrereqs, expName, runNumber, subRun, detectorName, dataLocation);
        vector<double> currentMonitorCounts;
        vector<double> currentGoodMacroCounts;

        for(auto& p : allCSPrereqs)
        {
            currentMonitorCounts.push_back(p.monitorCounts);
            currentGoodMacroCounts.push_back(p.goodMacroNumber);
        }

        runningMonitorCounts.push_back(currentMonitorCounts);
        runningGoodMacroCounts.push_back(currentGoodMacroCounts);
    }

    string outFileName = dataLocation + "/CSPrereqs.root";
    TFile* outFile = new TFile(outFileName.c_str(), "UPDATE");

    vector<TH1D*> runningFluxHistos;
    for(string targetName : config.target.TARGET_ORDER)
    {
        string runningFluxHistoName = targetName + "RunningFluxAvg";
        runningFluxHistos.push_back(new TH1D(runningFluxHistoName.c_str(),
                    runningFluxHistoName.c_str(), highSubrun-lowSubrun, lowSubrun, highSubrun));
    }

    for(int i=0; i<runningMonitorCounts.size(); i++)
    {
        for(int j=0; j<runningMonitorCounts[i].size(); j++)
        {
            double targetFluxAvg = runningMonitorCounts[i][j]/runningGoodMacroCounts[i][j];
            double blankFluxAvg = runningMonitorCounts[i][1]/runningGoodMacroCounts[i][1];

            runningFluxHistos[j]->SetBinContent(i+1, targetFluxAvg/blankFluxAvg);
        }
    }

    for(auto& histo : runningFluxHistos)
    {
        histo->Write();
    }

    outFile->Close();
}
