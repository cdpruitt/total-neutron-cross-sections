#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TMath.h"
#include "TLatex.h"
#include "TH1.h"
#include "TH2.h"

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

const int MAX_SUBRUN_NUMBER = 50;

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

    // Ingest data from every run in the run list
    int runNumber;

    vector<vector<CrossSection>> allCrossSections;

    vector<vector<CSPrereqs>> allCSPrereqs;

    int counter = 0;

    string line;
    while (runList >> line)
    {
        runNumber = stoi(line);

        // read in run config file
        config = Config(expName, runNumber);

        // Loop through all subruns of this run
        for(int subRun=0; subRun<=MAX_SUBRUN_NUMBER; subRun++)
        {
            vector<CSPrereqs> subRunCSPrereqs;

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

                cout << "Read subrun " << runNumber << " " << subRun << endl;

                subRunCSPrereqs.push_back(subRunData);
            }

            allCSPrereqs.push_back(subRunCSPrereqs);

            cout << "counter = " << counter++ << endl;

            CSPrereqs blank;
            for(auto& p : subRunCSPrereqs)
            {
                if(p.target.getName()=="blank" || p.target.getName()=="blankW")
                {
                    blank = p;
                }
            }

            if(blank.monitorCounts==0)
            {
                cerr << "Blank monitor counts was 0; skipping sub-run" << endl;
                continue;
            }

            vector<CrossSection> crossSections;
            for(auto& p : subRunCSPrereqs)
            {
                CrossSection cs;
                cs.calculateCS(p,blank);

                string blankCSFileName = "../" + dataLocation + "/literature.root"; 

                correctForBlank(cs, p, expName);
                crossSections.push_back(cs);;
            }

            string outFileName = dataLocation + "/subRunCS/" + to_string(runNumber) + "_" + to_string(subRun) + ".root";
            TFile* outFile = new TFile(outFileName.c_str(), "RECREATE");

            for(auto& cs : crossSections)
            {
                cs.createGraph(cs.name, cs.name);
            }

            allCrossSections.push_back(crossSections);

            outFile->Close();
        }
    }

    string outFileName = dataLocation + "/allSubRuns.root";

    TFile* outFile = new TFile(outFileName.c_str(), "RECREATE");

    vector<TMultiGraph*> multigraphs;

    for(auto& crossSection : allCrossSections[0])
    {
        multigraphs.push_back(new TMultiGraph());
    }

    for(int j=0; j<allCrossSections.size(); j++)
    {
        vector<CrossSection> crossSections = allCrossSections[j];

        for(int i=0; i<crossSections.size(); i++)
        {
            CrossSection cs = crossSections[i];

            TGraphErrors* t = new TGraphErrors(cs.getNumberOfPoints(),
                                      &cs.getEnergyValues()[0],
                                      &cs.getCrossSectionValues()[0],
                                      &cs.getEnergyErrors()[0],
                                      &cs.getCrossSectionErrors()[0]);

            string name = to_string(j) + "-" + to_string(i);
            t->SetNameTitle(name.c_str(),name.c_str());

            multigraphs[i]->Add(t);
        }
    }

    for(auto& multigraph : multigraphs)
    {
        multigraph->Write();
    }

    int numberOfCSPrereqs = 0;

    for(auto& subrun : allCSPrereqs)
    {
        for(auto& CSPrereqs : subrun)
        {
            numberOfCSPrereqs++;
        }
    }

    vector<TH1D*> monitorToDetectorRatios;

    for(int i=0; i<config.target.TARGET_ORDER.size(); i++)
    {
        string name = config.target.TARGET_ORDER[i] + "_mon/det";

        monitorToDetectorRatios.push_back(new TH1D(name.c_str(), name.c_str(), 
            numberOfCSPrereqs, 0, numberOfCSPrereqs));
    }

    int currentBin = 0;
    
    for(int i=0; i<allCSPrereqs.size(); i++)
    {
        vector<CSPrereqs> subRun = allCSPrereqs[i];

        CSPrereqs blank;
            
        for(auto csp : subRun)
        {
            if(csp.target.getName() == "blank")
            {
                blank = csp;
                break;
            }
        }

        for(int k=0; k<subRun.size(); k++)
        {
            CSPrereqs csp = subRun[k];

            for(int j=0; j<config.target.TARGET_ORDER.size(); j++)
            {
                if(csp.target.getName() == config.target.TARGET_ORDER[j])
                {
                    double ratio = csp.monitorCounts/csp.totalEventNumber;
                    monitorToDetectorRatios[j]->SetBinContent(currentBin, ratio);
                    break;
                }
            }
            currentBin++;

            for(int j=0; j<config.target.TARGET_ORDER.size(); j++)
            {
                if(csp.target.getName() == config.target.TARGET_ORDER[j])
                {
                    double ratio = csp.monitorCounts/csp.totalEventNumber;
                    monitorToDetectorRatios[j]->SetBinContent(currentBin, ratio);
                    break;
                }
            }
        }
    }

    for(auto& histo : monitorToDetectorRatios)
    {
        histo->Write();
    }

    /*vector<TH2D*> monitorToDetectorRatios2D;

    for(int i=0; i<config.target.TARGET_ORDER.size(); i++)
    {
        string name = config.target.TARGET_ORDER[i] + "_mon/det_2D";

        monitorToDetectorRatios2D.push_back(new TH2D(name.c_str(), name.c_str(),
            1000, 0, 50,
            1000, 0, 1));
    }

    for(int i=0; i<allCSPrereqs.size(); i++)
    {
        for(int j=0; j<config.target.TARGET_ORDER.size(); j++)
        {
            if(allCSPrereqs[i].target.getName() == config.target.TARGET_ORDER[j])
            {
                double countRatio = allCSPrereqs[i].totalEventNumber/allCSPrereqs[1].totalEventNumber;
                double fluxPerMacro = allCSPrereqs[i].monitorCounts/allCSPrereqs[i].goodMacroNumber;

                monitorToDetectorRatios2D[j]->Fill(fluxPerMacro,countRatio);
                break;
            }
        }
    }

    for(auto& histo : monitorToDetectorRatios2D)
    {
        histo->Write();
    }
    */

    vector<TH2D*> detectorToMonitor2Ds;

    for(int i=0; i<config.target.TARGET_ORDER.size(); i++)
    {
        string name = config.target.TARGET_ORDER[i] + "_det/mon_2D";

        detectorToMonitor2Ds.push_back(new TH2D(name.c_str(), name.c_str(),
            1000, 0.95, 1.05,
            1000, 0, 2));

        detectorToMonitor2Ds.back()->SetMarkerStyle(20);
    }

    for(int i=0; i<allCSPrereqs.size(); i++)
    {
        vector<CSPrereqs> subRun = allCSPrereqs[i];
        for(int k=0; k<subRun.size(); k++)
        {
            CSPrereqs csp = subRun[k];

            CSPrereqs blank;

            for(auto csp : subRun)
            {
                if(csp.target.getName() == "blank")
                {
                    blank = csp;
                    break;
                }
            }

            for(int j=0; j<config.target.TARGET_ORDER.size(); j++)
            {
                if(csp.target.getName() == config.target.TARGET_ORDER[j])
                {
                    double countRatio = csp.totalEventNumber/blank.totalEventNumber;
                    double macroRatio = (csp.monitorCounts/blank.monitorCounts)/(csp.goodMacroNumber/blank.goodMacroNumber);

                    detectorToMonitor2Ds[j]->Fill(macroRatio,countRatio);
                    break;
                }
            }
        }
    }

    for(auto& histo : detectorToMonitor2Ds)
    {
        histo->Write();
    }

    outFile->Close();
}
