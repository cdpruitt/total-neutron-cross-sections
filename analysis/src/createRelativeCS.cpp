/* How to use:
 *
 * ./createRelativeCS [first CS file name] [first CS graph name]
 *                    [second CS file name] [second CS graph name]
 *                    [output file name] [output graph name]
 */

#include "../include/dataSet.h"
#include "../include/crossSection.h"
#include "../include/plots.h"

#include "TFile.h"
#include "TGraphErrors.h"

#include <string>
#include <iostream>

using namespace std;

int main(int, char* argv[])
{
    string firstCSFileName = argv[1];
    string firstCSGraphName = argv[2];

    string secondCSFileName = argv[3];
    string secondCSGraphName = argv[4];

    string outputCSFileName = argv[5];
    string outputCSGraphName = argv[6];

    TFile* firstCSFile = new TFile(firstCSFileName.c_str(),"READ");
    TGraphErrors* firstCSGraph = (TGraphErrors*)firstCSFile->Get(firstCSGraphName.c_str());

    if(!firstCSGraph)
    {
        cerr << "Error: failed to find " << firstCSGraphName << " in " << firstCSGraphName << endl;
        exit(1);
    }

    TFile* secondCSFile = new TFile(secondCSFileName.c_str(),"READ");
    TGraphErrors* secondCSGraph = (TGraphErrors*)secondCSFile->Get(secondCSGraphName.c_str());

    if(!secondCSGraph)
    {
        cerr << "Error: failed to find " << secondCSGraphName << " in " << secondCSGraphName << endl;
        exit(1);
    }

    DataSet firstDS = DataSet(firstCSGraph,firstCSGraphName);
    DataSet secondDS = DataSet();

    for(int i=0; i<firstDS.getNumberOfPoints(); i++)
    {
        secondDS.addPoint(
                DataPoint(firstDS.getPoint(i).getXValue(),
                    firstDS.getPoint(i).getXError(),
                    secondCSGraph->Eval(firstDS.getPoint(i).getXValue()),
                    secondCSGraph->GetErrorY(firstDS.getPoint(i).getXValue()))); 
    }

    CrossSection firstCS = CrossSection();
    firstCS.addDataSet(firstDS);

    CrossSection secondCS = CrossSection();
    secondCS.addDataSet(secondDS);

    CrossSection relative = calculateRelative(firstCS, secondCS);

    TFile* outputCSFile = new TFile(outputCSFileName.c_str(),"UPDATE");
    relative.createCSGraph(outputCSGraphName.c_str(), outputCSGraphName.c_str());
    outputCSFile->Close();

    return 0;
}
