// this is a short 'script' to read cross section data from a ROOT graph
// and create a text version of the data in the graph

// how to use:
// ./readGraphToText path/to/ROOTfile graphName name/of/outputFile
//
// Example:
// ./readGraphToText /data2/analysis/total.root O18 O18.total

#include <iostream>
#include <fstream>
#include <sstream>

// ROOT classes
#include "TFile.h"

// Self-defined classes
#include "../include/dataSet.h"

using namespace std;

int main(int, char* argv[])
{
    string inputFileName = argv[1];
    string graphName = argv[2];
    string outputFileName = argv[3];

    TFile* file = new TFile(inputFileName.c_str(), "READ");

    if(!file)
    {
        cerr << "Error: could not open file " << inputFileName << endl;
        return 1;
    }

    TGraphAsymmErrors* graph = (TGraphAsymmErrors*)file->Get(graphName.c_str());
    if(!graph)
    {
        cerr << "Error: could not find graph " << graphName << " in file " << inputFileName
            << endl;
        return 1;
    }

    DataSet graphData(graph, graphName);

    ofstream outputFile(outputFileName.c_str());

    outputFile << graphData.getReference() << endl;
    outputFile << "MeV        barns        error(barns)" << endl;

    for(int i=0; i<graphData.getNumberOfPoints(); i++)
    {
        DataPoint p = graphData.getPoint(i);

        outputFile << p.getXValue() << "        " << p.getYValue()
            << "        " << p.getYError() << endl;
    }

    file->Close();

    return 0;
}
