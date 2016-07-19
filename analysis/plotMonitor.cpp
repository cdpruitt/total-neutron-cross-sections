#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMath.h"

using namespace std;

const int NUMBER_OF_TARGETS = 6;
const int MAX_SUBRUN_NUMBER = 30;

const vector<string> targetNames = {"blankDPP", "shortCarbonDPP", "longCarbonDPP", "Sn112DPP", "NatSnDPP", "Sn124DPP"}; 
const vector<string> targetNamesWaveform = {"blankWaveform", "shortCarbonWaveform", "longCarbonWaveform", "Sn112Waveform", "NatSnWaveform", "Sn124Waveform"}; 

int main(int argc, char *argv[])
{
    string runNumber = argv[1];
    string driveName = argv[2];

    // create vectors for holding cross section data:
    // crossSections[target number]->at(subrun number)->at(data point)
    vector<vector<double>*> monitorCounts;
    vector<double> subrunNumber;

    ifstream monitorCountFile;

    // Loop over subruns to create average histos
    for(int i = 0; i<=MAX_SUBRUN_NUMBER; i++)
    {
        stringstream run;
        run << i;
        string runString = run.str();

        // Calculate the name of the next sub-run
        if(i < 10)
        {
            monitorCountFile.open(driveName + "/analysis/run" + runNumber + "/run" + runNumber + "-000" + runString + "_monitorCounts.log");
        }

        else if(i < 100)
        {
            monitorCountFile.open(driveName + "/analysis/run" + runNumber + "/run" + runNumber + "-00" + runString + "_monitorCounts.log");
        }

        else
        {
            cout << "Error: subrun number too large." << endl;
            exit(1);
        }

        if(!monitorCountFile.is_open())
        {
            cout << "Can't open subrun " << runNumber << " " << i << endl;
            continue;
        }

        // Successfully found the sub-run of interest
        cout << "Adding run " << runNumber << " " << i <<endl;

        subrunNumber.push_back(i);
        monitorCounts.push_back(new vector<double>);
        
        // Pull out the monitor data
        while(!monitorCountFile.eof())
        {
            double temp;
            monitorCountFile >> temp;
            monitorCounts.back()->push_back(temp);
        }

        monitorCountFile.close();
        // End of loop - move to next sub-run
    }

    // Create output file to contain summed histos
    stringstream outfileName;
    outfileName << driveName << "/analysis/run" << runNumber << "/" << "monitor.root";
    TFile *outfile = new TFile(outfileName.str().c_str(),"RECREATE");

    vector<vector<double>*> monitorRemapped;
    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        monitorRemapped.push_back(new vector<double>);
        for(int j=0; j<subrunNumber.size(); j++)
        {
            monitorRemapped[i]->push_back(monitorCounts[j]->at(i));
            if(i==0)
            {
                cout << monitorRemapped[i]->at(j) << endl;
            }
        }
    }

    vector<vector<double>*> monitorSums;
    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        monitorSums.push_back(new vector<double>);
        monitorSums.back()->resize(subrunNumber.size());
        for(int j=0; j<subrunNumber.size(); j++)
        {
            for(int k=0; k<=j; k++)
            {
                    monitorSums[i]->at(j) += monitorRemapped[i]->at(k);
            }
        }
    }

    vector<vector<double>*> monitorRatio;
    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        monitorRatio.push_back(new vector<double>);
        monitorRatio.back()->resize(subrunNumber.size());
        for(int j=0; j<subrunNumber.size(); j++)
        {
            if(monitorSums[0]->at(j)!=0)
            {
                monitorRatio[i]->at(j) = monitorSums[i]->at(j)/monitorSums[0]->at(j);
            }
            if(i==0)
            {
            //    cout << monitorRatio[i]->at(j) << endl;
            }
        }
    }

    // create new graphs to display the monitor counts over time
    for(int i=0; i<monitorRemapped.size(); i++)
    {
        TGraph* graph = new TGraph(subrunNumber.size(),&subrunNumber[0],&monitorRatio[i]->at(0));
        graph->SetNameTitle(targetNames[i].c_str(),targetNames[i].c_str());
        graph->Write();
    }

    outfile->Write();
    outfile->Close();
}
