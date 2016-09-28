// this is a short 'script' to read cross section data from a text file
// and create ROOT graphs of the data
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
    string outFileName = inputDirectory + outputFile;

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
        dummy = inputDirectory + dummy;
        fileNames.push_back(dummy);
    }

    // open/create output file
    TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

    vector<DataSet*> allData;

    for(string s : fileNames)
    {
        allData.push_back(new DataSet(s));
    }

    /*if(allData.size()>1)
    {
        DataSet Sum124_112 = allData[2].plus(allData[3],"Sn124+Sn112");
        DataSet Diff124_112 = allData[2].minus(allData[3],"Sn124-Sn112");
        DataSet RelDiff124_112 = Diff124_112.divideBy(Sum124_112,"(Sn124-Sn112)/(Sn124+Sn112)");
    }*/

    // clean up
    //outFile->Write();
    outFile->Close();

    return 0;
}
