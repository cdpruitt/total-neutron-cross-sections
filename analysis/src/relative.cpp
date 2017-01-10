// how to use:
// ./relative experimentName path/to/literature/data.root path/to/relative/crossSections.root

#include <iostream>
#include <fstream>
#include <string>
#include "TFile.h"
#include "TGraphErrors.h"

#include "../include/experiment.h"
#include "../include/dataSet.h"
#include "../include/crossSection.h"

using namespace std;

int main(int, char* argv[])
{
    string expName = argv[1];
    string litDataLocation = argv[2];
    string relativeFileLocation = argv[3];

    vector<pair<string,string>> relativeTargets = getRelativePlotNames(expName,"litRelativePlots.txt");

    ifstream litData(litDataLocation.c_str());
    if(!litData.is_open())
    {
        cerr << "Attempted to open literature data, but failed to find " << litDataLocation << endl;
        exit(1);
    }
    litData.close();

    TFile* litDataFile = new TFile(litDataLocation.c_str(),"READ");
    TFile* relativeFile = new TFile(relativeFileLocation.c_str(),"UPDATE");

    for(pair<string,string> p : relativeTargets)
    {
        cout << "Producing relative cross section plot of " << p.first << " and " << p.second << endl;

        TGraphErrors* firstGraph = (TGraphErrors*)litDataFile->Get(p.first.c_str());
        if(!firstGraph)
        {
            cerr << "Error: failed to open " << p.first << " graph for relative cross section calculation." << endl;
            exit(1);
        }

        TGraphErrors* secondGraph = (TGraphErrors*)litDataFile->Get(p.second.c_str());
        if(!secondGraph)
        {
            cerr << "Error: failed to open " << p.second << " graph for relative cross section calculation." << endl;
            exit(1);
        }

        CrossSection first;
        first.addDataSet(DataSet(firstGraph,p.first));

        CrossSection second;
        DataSet secondDataSet;

        // for each y-value of the first CS graph, read the y-value of the
        // second graph and the y-error

        DataSet firstDataSet = first.getDataSet();
        for(int i=0; i<firstDataSet.getNumberOfPoints(); i++)
        {
            secondDataSet.addPoint(
                    DataPoint(firstDataSet.getPoint(i).getXValue(),
                        firstDataSet.getPoint(i).getXError(),
                        secondGraph->Eval(firstDataSet.getPoint(i).getXValue()),
                        secondGraph->GetErrorY(firstDataSet.getPoint(i).getXValue()))); 
        }

        second.addDataSet(secondDataSet);

        CrossSection relative = calculateRelative(first,second);
        string relativeName = "#frac{#sigma_{" + p.first +
                              "}-#sigma_{" + p.second +
                              "}}{#sigma_{" + p.first +
                              "}+#sigma_{" + p.second + "}}";
        relative.createCSGraph(relativeName.c_str());
        //cout << "Relative plot " << relative.getDataSet().getReference() <<
        //            " RMS error: " << relative.calculateRMSError() << endl;
    }

    litDataFile->Close();
    relativeFile->Close();

    return 0;
}
