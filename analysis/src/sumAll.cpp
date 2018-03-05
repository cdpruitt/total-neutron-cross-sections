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
#include "../include/correctForBackground.h"

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
    vector<vector<double>> runningFluxAvg;

    bool printRunCount;

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
            printRunCount = false;

            // Loop through all target positions in this subrun
            for(int j=0; (size_t)j<config.target.TARGET_ORDER.size(); j++)
            {
                // pull data needed for CS calculation from subrun 
                string targetDataLocation = "../" + expName + "/targetData/" + config.target.TARGET_ORDER[j] + ".txt";
                CSPrereqs subRunData(targetDataLocation);

                if(readSubRun(subRunData, expName, runNumber, subRun, detectorName, dataLocation))
                {
                    break;
                }

                printRunCount = true;

                // find the correct CSPrereqs to add this target's data to

                for(CSPrereqs& csp : allCSPrereqs)
                {
                    if(csp.target.getName() == subRunData.target.getName())
                    {
                        // add subrun data to total
                        csp = csp + subRunData;
                    }
                }

                vector<double> currentFluxAvg;

                for(auto& p : allCSPrereqs)
                {
                    currentFluxAvg.push_back(p.monitorCounts);
                }

                runningFluxAvg.push_back(currentFluxAvg);

            }

            if(printRunCount)
            {
                cout << "Read " << runNumber << "-" << subRun << endl;
            }
        }
    }

    string outFileName = dataLocation + "/total.root";
    TFile* outFile = new TFile(outFileName.c_str(), "UPDATE");

    CSPrereqs blank;
    for(auto& p : allCSPrereqs)
    {
        correctForBackground(p);

        string energyHistoName = p.target.getName();
        energyHistoName = energyHistoName + "Energy";

        p.energyHisto = convertTOFtoEnergy(p.TOFHisto, energyHistoName.c_str());

        if(p.target.getName()=="blank" || p.target.getName()=="blankW")
        {
            blank = p;
        }

        cout << "all CS Prereqs name = " << p.target.getName() << endl;

        double totalCounts = 0;
        int numberOfBins = p.TOFHisto->GetNbinsX();

        for(int i=1; i<=numberOfBins; i++)
        {
            double tempCounts = p.TOFHisto->GetBinContent(i);
            if(tempCounts < 0)
            {
                continue;
            }

            totalCounts += tempCounts;
        }

        cout << endl << "Total statistics over all runs: " << endl << endl;
        cout << p.target.getName() << ": total events in TOF histo = "
            << totalCounts << ", total monitor events = "
            << p.monitorCounts << ", good macro number = "
            << p.goodMacroNumber << ", total macro number = "
            << p.totalMacroNumber << ", total event number = "
            << p.totalEventNumber << endl;

        string name = p.target.getName() + "TOF";
        p.TOFHisto->SetNameTitle(name.c_str(),name.c_str());
        p.TOFHisto->SetDirectory(outFile);
        p.TOFHisto->Write();

        name = p.target.getName() + "TOFUncorrected";
        p.uncorrectedTOFHisto->SetNameTitle(name.c_str(),name.c_str());
        p.uncorrectedTOFHisto->SetDirectory(outFile);
        p.uncorrectedTOFHisto->Write();
    }

    vector<CrossSection> crossSections;
    for(auto& p : allCSPrereqs)
    {
        CrossSection cs;
        cs.calculateCS(p,blank);

        correctForBlank(cs, p, expName);
        crossSections.push_back(cs);;
    }

    outFile->cd();

    for(auto& cs : crossSections)
    {
        cs.createGraph(cs.name, cs.name);
    }

    outFile->Close();

    // read literature data and bin to appropriate energy range
    string litDirectory = "../" + expName + "/literatureData";
    string litOutputName = dataLocation + "/literatureData.root";
    if(readLitData(litDirectory, litOutputName, config))
    {
        cerr << "Error: failed to produce properly binned literature cross sections. Exiting..." << endl;
        return 1;
    }

    cout << "Finished binning/plotting literature cross sections." << endl;

    return 0;
}
