/* How to use:
 *
 * ./produceRunningRMS [exp CS file name] [exp CS graph name]
 *                     [lit CS file name] [lit CS graph name]
 *                     [output file name] [output graph name]
 */

#include "../include/dataSet.h"
#include "../include/crossSection.h"
#include "../include/plots.h"
#include "../include/config.h"
#include "../include/CSUtilities.h"

#include "TFile.h"
#include "TGraphAsymmErrors.h"

#include <string>
#include <iostream>

using namespace std;

Config config;

int main(int, char* argv[])
{
    string expCSFileName = argv[1];
    string expCSGraphName = argv[2];
    string litCSFileName = argv[3];
    string litCSGraphName = argv[4];
    string outputFileName = argv[5];
    string outputGraphName = argv[6];

    TFile* expCSFile = new TFile(expCSFileName.c_str(),"READ");
    TGraphAsymmErrors* expCSGraph = (TGraphAsymmErrors*)expCSFile->Get(expCSGraphName.c_str());
    if(!expCSGraph)
    {
        cerr << "Error: failed to find " << expCSGraphName << " in " << expCSFileName << endl;
        exit(1);
    }

    TFile* litCSFile = new TFile(litCSFileName.c_str(),"READ");
    TGraphAsymmErrors* litCSGraph = (TGraphAsymmErrors*)litCSFile->Get(litCSGraphName.c_str());
    if(!litCSGraph)
    {
        cerr << "Error: failed to find " << litCSGraphName << " in " << litCSFileName << endl;
        exit(1);
    }

    DataSet expDataSet = DataSet(expCSGraph, expCSGraphName);
    DataSet litDataSet = DataSet();

    // for each value of the expCS graph, read the y-value of the litCS graph
    // and the y-error
    for(int i=0; i<expDataSet.getNumberOfPoints(); i++)
    {
        litDataSet.addPoint(
                DataPoint(expDataSet.getPoint(i).getXValue(),
                    expDataSet.getPoint(i).getXError(),
                    litCSGraph->Eval(expDataSet.getPoint(i).getXValue()),
                    litCSGraph->GetErrorY(i))); 
    }

    expCSFile->Close();
    litCSFile->Close();

    TFile* outputFile = new TFile(outputFileName.c_str(),"UPDATE");

    produceRunningRMS(expDataSet, litDataSet, outputGraphName);

    return 0;
}
