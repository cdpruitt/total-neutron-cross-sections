// this is a short 'script' to read cross section data from a text file
// and create ROOT graphs of the data

// how to use:
// ./readLitData path/to/literatureData name/of/outputFile

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
    string inputDirectory = argv[1];
    string outputFile = argv[2];
    string inFileName = inputDirectory + "/filesToRead.txt";
    string outFileName = inputDirectory + "/" + outputFile;

    ifstream inFile(inFileName);
    if(!inFile.is_open())
    {
        cout << "Failed to open " << inFileName << endl;
        exit(1);
    }

    string dummy;
    vector<string> fileNames;

    while(inFile >> dummy)
    {
        dummy = inputDirectory + "/" + dummy;
        fileNames.push_back(dummy);
    }

    // recreate output file
    TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

    vector<DataSet> allData;

    for(string s : fileNames)
    {
        cout << "Creating plot for " << s << endl;
        allData.push_back(DataSet(s));
        outFile->cd();
        allData.back().getPlot()->Write();
    }

    outFile->Close();

    return 0;
}
