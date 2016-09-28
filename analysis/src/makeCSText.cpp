// this is a short 'script' to turn a ROOT graph of cross sections into a
// text file, formatted for inclusion in DOM fits
#include <iostream>
#include "TFile.h"
#include <string>
#include <include/dataPoint.h>

void makeCSText(string inputFileName, string graphName, string outputFileName)
{
    // open cross section file
    TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
    if(!inputFile)
    {
        cout << "Failed to find cross section file in makeCSText. Exiting..." << endl;
        exit(1);
    }

    TGraphErrors* plot = (TGraphErrors*)gDirectory->Get(graphName.c_str());
    if(!plot)
    {
        cout << "Failed to open cross section graph in makeCSText. Exiting..." << endl;
        exit(1);
    }

    int n = plot->GetN();
    double x;
    double y;
    double yError;

    // create output file
    ofstream outputFile(outputFileName);

    for(int i=0; i<n; i++)
    {
        plot->GetPoint(i,x,y);
        yError = plot->GetErrorY(i);
        outputFile << x << "    " << y << "    " << yError << endl;
    }

    inputFile->Close();
    outputFile.close();
}
