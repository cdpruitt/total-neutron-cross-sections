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
    
    // read in run config file
    config = Config(expName, runNumber);

    readTargetData(allCSPrereqs, expName);

    // Loop through all subruns of this run
    for(int subRun=lowSubrun; subRun<=highSubrun; subRun++)
    {
        readSubRun(allCSPrereqs, expName, runNumber, subRun, detectorName, dataLocation);
    }

    string outFileName = dataLocation + "/total.root";
    TFile* outFile = new TFile(outFileName.c_str(), "UPDATE");

    CSPrereqs blank;
    for(auto& p : allCSPrereqs)
    {
        if(p.target.getName()=="blank")
        {
            blank = p;
            break;
        }
    }

    if(blank.target.getName()!="blank")
    {
        cerr << "Error: failed to find blank target for cross section calculation." << endl;
        return 1;
    }

    for(auto& p : allCSPrereqs)
    {
        long totalCounts = 0;
        for(int i=0; i<p.energyHisto->GetNbinsX(); i++)
        {
            long tempCounts = p.energyHisto->GetBinContent(i);
            if(tempCounts < 0)
            {
                continue;
            }

            totalCounts += tempCounts;
        }

        cout << endl << "Total statistics over all runs: " << endl << endl;
        cout << p.target.getName() << ": total events in energy histo = "
            << totalCounts << ", total monitor events = "
            << p.monitorCounts << ", good macro number = "
            << p.goodMacroNumber << ", total macro number = "
            << p.totalMacroNumber << endl;

        p.energyHisto->SetDirectory(outFile);
        p.energyHisto->Write();

        string name = p.target.getName() + "TOF";
        p.TOFHisto->SetNameTitle(name.c_str(),name.c_str());
        p.TOFHisto->SetDirectory(outFile);

        p.TOFHisto->Write();
    }

    vector<CrossSection> crossSections;
    for(auto& p : allCSPrereqs)
    {
        crossSections.push_back(CrossSection());
        crossSections.back().calculateCS(p,blank);

        cout << "Created cross section for " << crossSections.back().name << " target." << endl;
    }

    for(auto& cs : crossSections)
    {
        cs.createGraph(cs.name, cs.name);
    }

    outFile->Close();
}
