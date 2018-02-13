/* How to use:
 *
 * ./mergeCS [first CS filename] [first CS graph name]
 *           [second CS filename] [second CS graph name]
 *           [juncture point] [output CS filename] [outputCS graph name]
 */

#include "../include/dataSet.h"
#include "../include/crossSection.h"
#include "../include/plots.h"
#include "../include/config.h"
#include "../include/CSUtilities.h"

#include "TFile.h"
#include "TGraphErrors.h"

#include <string>
#include <iostream>

using namespace std;

Config config;

int main(int, char* argv[])
{
    string firstCSFileName = argv[1];
    string firstCSGraphName = argv[2];
    string secondCSFileName = argv[3];
    string secondCSGraphName = argv[4];
    double juncture = stod(argv[5]);
    string outputCSFileName = argv[6];
    string outputCSGraphName = argv[7];

    TFile* firstCSFile = new TFile(firstCSFileName.c_str(),"READ");
    if(!firstCSFile->IsOpen())
    {
        cerr << "Error: couldn't open " << firstCSFileName << "." << endl;
        return 1;
    }

    TFile* secondCSFile = new TFile(secondCSFileName.c_str(),"READ");
    if(!secondCSFile->IsOpen())
    {
        cerr << "Error: couldn't open " << secondCSFileName << "." << endl;
        firstCSFile->Close();
        return 1;
    }

    TGraphErrors* firstCSGraph = (TGraphErrors*)firstCSFile->Get(firstCSGraphName.c_str());
    if(!firstCSGraph)
    {
        cerr << "Error: couldn't find " << firstCSGraphName << " in " << firstCSFileName << "." << endl;
        firstCSFile->Close();
        secondCSFile->Close();
        return 1;
    }

    TGraphErrors* secondCSGraph = (TGraphErrors*)secondCSFile->Get(secondCSGraphName.c_str());
    if(!secondCSGraph)
    {
        cerr << "Error: couldn't find " << secondCSGraphName << " in " << secondCSFileName << "." << endl;
        firstCSFile->Close();
        secondCSFile->Close();
        return 1;
    }

    CrossSection firstCS;
    firstCS.addDataSet(DataSet(firstCSGraph,firstCSGraph->GetName()));

    CrossSection secondCS;
    secondCS.addDataSet(DataSet(secondCSGraph,secondCSGraph->GetName()));

    CrossSection mergedCS = mergeCrossSections(firstCS, juncture, secondCS);

    firstCSFile->Close();
    secondCSFile->Close();

    TFile* outputCSFile = new TFile(outputCSFileName.c_str(),"UPDATE");
    mergedCS.createGraph(outputCSGraphName, outputCSGraphName);

    outputCSFile->Close();

    return 0;
}
