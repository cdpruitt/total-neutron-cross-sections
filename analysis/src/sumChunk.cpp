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
#include "../include/correctForBackground.h"

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

            // find the correct CSPrereqs to add this target's data to

            for(CSPrereqs& csp: allCSPrereqs)
            {
                if(csp.target.getName() == subRunData.target.getName())
                {
                    // add subrun data to total
                    csp = csp + subRunData;
                }
            }
        }
    }

    string outFileName = dataLocation + "/total.root";
    TFile* outFile = new TFile(outFileName.c_str(), "UPDATE");

    CSPrereqs blank;
    for(auto& p : allCSPrereqs)
    {
        //correctForBackground(p);

        string energyHistoName = p.target.getName();
        energyHistoName = energyHistoName + "Energy";

        p.energyHisto = convertTOFtoEnergy(p.TOFHisto, energyHistoName.c_str());

        if(p.target.getName()=="blank" || p.target.getName()=="blankW")
        {
            blank = p;
        }

        double totalCounts = 0;
        int numberOfBins = p.energyHisto->GetNbinsX();

        for(int i=1; i<=numberOfBins; i++)
        {
            double tempCounts = p.energyHisto->GetBinContent(i);
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
            << p.totalMacroNumber << ", total event number = "
            << p.totalEventNumber << endl;

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
        CrossSection cs = CrossSection();
        cs.calculateCS(p,blank);
        correctForBlank(cs, p, expName);

        crossSections.push_back(cs);

        cout << "Created cross section for " << crossSections.back().name << " target." << endl;
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
