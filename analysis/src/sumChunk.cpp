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

    producePlots(dataLocation, allCSPrereqs);
}
