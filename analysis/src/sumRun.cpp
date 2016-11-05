#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TMath.h"

#include "../include/dataSet.h"

using namespace std;

const int MAX_SUBRUN_NUMBER = 99;

//const vector<string> targetNames = {"shortCarbon", "longCarbon", "Sn112", "NatSn", "Sn124"}; 
const vector<string> targetNames = {"NatPb", "longCarbon", "Sn112", "NatSn", "Sn124"};

// Extract each point from a graph and store their positions in two vectors,
// xValues and yValues
DataSet extractGraphData(TGraphErrors* graph, string name)
{
    int numPoints = graph->GetN();
    
    vector<double> xValues(numPoints);
    vector<double> yValues(numPoints);
    vector<double> yErrors(numPoints);

    for(int k=0; k<numPoints; k++)
    {
        graph->GetPoint(k,xValues[k],yValues[k]);
        yErrors[k] = graph->GetErrorY(k);
    }

    return DataSet(xValues,yValues,yErrors,name);
}

bool fileExists(string fileName)
{
    ifstream infile(fileName);
    return infile.good();
}

vector<DataSet> readGraphs(
        string runNumber,
        string driveName,
        string fileType,
        string target
        )
{

    // hold each graph's data for this target
    vector<DataSet> dataSets;

    TFile* infile;

    // Loop over subruns to read data 
    for(int i = 0; i<=MAX_SUBRUN_NUMBER; i++)
    {
        stringstream inFileName;

        // Calculate the name of the next sub-run
        if(i < 10)
        {
            inFileName << driveName << "/analysis/" << runNumber << "/000"
                       << i << "/" << fileType << ".root";
        }

        else if(i < 100)
        {
            inFileName << driveName << "/analysis/" << runNumber << "/00"
                       << i << "/" << fileType << ".root";
        }

        // Attempt to open the sub-run
        if(!fileExists(inFileName.str()))
        {
            continue;
        }

        infile = new TFile(inFileName.str().c_str());
        cout << "Adding target " << target << " in run " << runNumber << " " << i << endl;

        // Open the graph
        TGraphErrors * graph = (TGraphErrors*)infile->Get(target.c_str());
        if(!graph)
        {
            cout << "Couldn't find graph for " << target << endl;
            continue;
        }

        dataSets.push_back(extractGraphData(graph,target));

        infile->Close();
        // End of loop - move to next sub-run
    }

    return dataSets;
}

void produceAverages(vector<DataSet> targetData)
{
    DataSet totalData = targetData[0];
    string name = totalData.getReference();

    for(int i=1; (size_t)i<targetData.size(); i++)
    {
        totalData = totalData + targetData[i];
    }

    DataSet averageData = totalData/targetData.size();
    TGraphErrors* plot = averageData.createPlot(name);
}

int main(int, char *argv[])
{
    string runNumber = argv[1];
    string driveName = argv[2];

    vector<vector<DataSet>> runData;

    for(string t : targetNames)
    {
        runData.push_back(readGraphs(runNumber, driveName, "cross-sections", t));
    }

    // Create output file to contain summed histos
    stringstream outfileName;
    outfileName << driveName << "/analysis/" << runNumber << "/" << "sum.root";
    TFile *outfile = new TFile(outfileName.str().c_str(),"RECREATE");

    gDirectory->mkdir("high_threshold","high_threshold");
    gDirectory->cd("high_threshold");

    for(vector<DataSet> target : runData)
    {
        produceAverages(target);
    }

    outfile->Write();
    outfile->Close();
}
