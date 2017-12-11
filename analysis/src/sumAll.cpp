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
#include "../include/crossSection.h"
#include "../include/experiment.h"
#include "../include/plots.h"
#include "../include/CSUtilities.h"

using namespace std;

Config config;

const int MAX_SUBRUN_NUMBER = 100;

int main(int, char* argv[])
{
    string dataLocation = argv[1];

    string expName = argv[2]; // experiment directory where runs to-be-sorted
    // are listed

    string detectorName = argv[3]; // detector name to be used for calculating cross sections

    // Open run list
    string runListName = "../" + expName + "/runsToSort.txt";
    ifstream runList(runListName);
    if(!runList.is_open())
    {
        cerr << "Error: couldn't find runlist at " << runListName << endl;
        exit(1);
    }
    cout << endl;

    // store run data in a "cross section prerequisites" structure
    vector<CSPrereqs> allCSPrereqs;

    // Ingest data from every run in the run list
    int runNumber;

    string line;
    while (runList >> line)
    {
        runNumber = stoi(line);

        // read in run config file
        config = Config(expName, runNumber);

        readTargetData(allCSPrereqs, expName);

        // Loop through all subruns of this run
        for(int subRun=0; subRun<=MAX_SUBRUN_NUMBER; subRun++)
        {
            readSubRun(allCSPrereqs, expName, runNumber, subRun, detectorName, dataLocation);
        }
    }

    string outFileName = dataLocation + "/total.root";
    TFile* outFile = new TFile(outFileName.c_str(), "UPDATE");

    CSPrereqs blank;
    for(auto& p : allCSPrereqs)
    {
        if(p.target.getName()=="blank")
        {
            blank = p;
        }

        cout << "all CS Prereqs name = " << p.target.getName() << endl;

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
        CrossSection cs;
        cs.calculateCS(p,blank);
        crossSections.push_back(cs);;
    }

    for(auto& cs : crossSections)
    {
        cs.createGraph(cs.name, cs.name);
    }

    outFile->Close();
}
