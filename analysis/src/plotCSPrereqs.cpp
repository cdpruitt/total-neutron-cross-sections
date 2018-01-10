#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
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
    vector<vector<double>> runningFluxAvg;
    
    // read in run config file
    config = Config(expName, runNumber);

    readTargetData(allCSPrereqs, expName);

    // Loop through all subruns of this run
    for(int subRun=lowSubrun; subRun<=highSubrun; subRun++)
    {
        readSubRun(allCSPrereqs, expName, runNumber, subRun, detectorName, dataLocation);

        vector<double> currentFluxAvg;

        for(auto& p : allCSPrereqs)
        {
            currentFluxAvg.push_back(p.monitorCounts);
        }

        runningFluxAvg.push_back(currentFluxAvg);
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

    for(unsigned int i=0; i<runningFluxAvg.size(); i++)
    {
        for(unsigned int j=0; j<runningFluxAvg[i].size(); j++)
        {
            runningFluxHistos[j]->SetBinContent(i+1, runningFluxAvg[i][1]/runningFluxAvg[i][j]);
        }
    }

    for(auto& histo : runningFluxHistos)
    {
        histo->Write();
    }

    outFile->Close();
}
