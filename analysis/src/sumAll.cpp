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

    producePlots(dataLocation, allCSPrereqs);
}
