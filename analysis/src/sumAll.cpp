#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TLatex.h"

#include "../include/target.h"
#include "../include/dataSet.h"
#include "../include/dataPoint.h"
#include "../include/CSPrereqs.h"
#include "../include/analysisConstants.h"
#include "../include/crossSection.h"

using namespace std;

const int MAX_SUBRUN_NUMBER = 200;

struct Plots
{
    vector<TGraphErrors*> CSGraphs;
    vector<TGraphErrors*> CSGraphsWaveform;

    vector<TGraphErrors*> relativeCSGraphs;
    vector<TGraphErrors*> relativeCSGraphsWaveform;
} plots;

// Extract each point from a graph and store their positions in two vectors,
// xValues and yValues
void extractGraphData(
        TGraphErrors* graph,
        vector<double>* xValues,
        vector<double>* xError,
        vector<double>* yValues,
        vector<double>* yError)
{
    int numPoints = graph->GetN();
    
    xValues->resize(numPoints);
    yValues->resize(numPoints);
    xError->resize(numPoints);
    yError->resize(numPoints);

    for(int k=0; k<numPoints; k++)
    {
        graph->GetPoint(k,xValues->at(k),yValues->at(k));
        xError->at(k) = graph->GetErrorX(k);
        yError->at(k) = graph->GetErrorY(k);
    }
}

// Extract each point from a graph and store their positions in a DataSet
void extractGraphData(
        TGraphErrors* graph,
        DataSet& dataSet)
{
    int numPoints = graph->GetN();

    double xValue;
    double xError;
    double yValue;
    double yError;

    for(int k=0; k<numPoints; k++)
    {
        graph->GetPoint(k,xValue,yValue);
        xError = graph->GetErrorX(k);
        yError = graph->GetErrorY(k);

        DataPoint dataPoint(xValue, xError, yValue, yError);
        dataSet.addPoint(dataPoint);
    }
}


// Calculate the root-mean-squared difference between two vectors
double calculateRMS(vector<double> graph1Data, vector<double> graph2Data)
{
    double rms = 0;
    for(int i=0; (size_t)i<graph1Data.size(); i++)
    {
        rms += pow(graph1Data.at(i)-graph2Data.at(i),2);
    }
    rms /= graph1Data.size();
    rms = pow(rms,0.5);
    return rms;
}

// Scale a dataset's y-values by the ratio of two other datasets
DataSet scale(DataSet setToScale, DataSet expReference, DataSet litReference)
{
    DataSet scaleFactor;

    int n = setToScale.getNumberOfPoints();
    for(int i=0; i<n; i++)
    {
        DataPoint p = litReference.getPoint(i)/expReference.getPoint(i);
        scaleFactor.addPoint(p);
    }

    DataSet scaledSet = scaleFactor*setToScale;

    return scaledSet;
}

vector<string> getTargetNames(string expName)
{
    string targetNamesLocation = "../" + expName + "/targetNames.txt";
    ifstream dataFile(targetNamesLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Failed to find target names in " << targetNamesLocation << std::endl;
        exit(1);
    }

    string str;
    vector<string> targetNames;

    while(getline(dataFile,str))
    {
        targetNames.push_back(str);
    }

    return targetNames;
}

// Determine the target order for a given run
vector<string> getTargetOrder(string expName, int runNumber)
{
    string targetOrderLocation = "../" + expName + "/targetOrder.txt";
    ifstream dataFile(targetOrderLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Failed to find target order data in " << targetOrderLocation << std::endl;
        exit(1);
    }

    string str;
    vector<string> targetOrder;

    while(getline(dataFile,str))
    {
        // ignore comments in data file
        string delimiter = "-";
        string token = str.substr(0,str.find(delimiter));
        if(!atoi(token.c_str()))
        {
            // This line starts with a non-integer and is thus a comment; ignore
            continue;
        }

        // parse data lines into space-delimited tokens
        vector<string> tokens;
        istringstream iss(str);
        copy(istream_iterator<string>(iss),
                istream_iterator<string>(),
                back_inserter(tokens));

        // extract run numbers from first token
        string lowRun = tokens[0].substr(0,tokens[0].find(delimiter));
        tokens[0] = tokens[0].erase(0,tokens[0].find(delimiter) + delimiter.length());

        delimiter = "\n";
        string highRun = tokens[0].substr(0,tokens[0].find(delimiter));
        
        if(atoi(lowRun.c_str()) <= runNumber && runNumber <= atoi(highRun.c_str()))
        {
            for(int i=1; (size_t)i<tokens.size(); i++)
            {
                targetOrder.push_back(tokens[i]);
            }
            break;
        }
    }

    return targetOrder;
}

int main(int, char* argv[])
{
    string dataLocation = argv[1]; // name of directory where analysis is stored
                                   // (omit trailing slash)
    string expName = argv[2];      // experiment directory where runs to-be-sorted
                                   // are listed
    string ROOTFileName = argv[3]; // name of ROOT files that contain data
                                   // used to calculate cross sections

    // Create a CSPrereqs for each target to hold data from all the runs
    vector<CSPrereqs> allData;

    vector<string> targetNames = getTargetNames(expName);
    for(string targetName : targetNames)
    {
        string targetDataLocation = "../" + expName + "/targetData/" + targetName + ".txt";
        allData.push_back(CSPrereqs(targetDataLocation));
    }

    // Open runlist
    string runListName = "../" + expName + "/runsToSort.txt";
    ifstream runList(runListName);
    if(!runList.is_open())
    {
        cerr << "Error: couldn't find runlist at " << runListName << endl;
        exit(1);
    }

    cout << endl;

    // Runlist open - loop through all runs
    string runNumber;
    while (runList >> runNumber)
    {
        // Loop through all subruns of this run
        for(int subRun=0; subRun<=MAX_SUBRUN_NUMBER; subRun++)
        {
            stringstream subRunFormatted;
            subRunFormatted << setfill('0') << setw(4) << subRun;

            // open subrun
            stringstream inFileName;
            inFileName << dataLocation << "/" << runNumber << "/"
                       << subRunFormatted.str() << "/" << ROOTFileName << ".root";
            ifstream f(inFileName.str());
            if(!f.good())
            {
                // failed to open this sub-run - skip to the next one
                continue;
            }

            f.close();

            TFile* inFile = new TFile(inFileName.str().c_str(),"READ");

            // get target order for this run
            vector<string> targetOrder = getTargetOrder(expName, stoi(runNumber));

            cout << "Adding " << runNumber << ", subrun " << subRun << "\r";
            fflush(stdout);

            // Loop through all target positions in this subrun
            for(int j=0; (size_t)j<targetOrder.size(); j++)
            {
                // pull data needed for CS calculation from subrun 
                string targetDataLocation = "../" + expName + "/targetData/" + targetOrder[j] + ".txt";
                CSPrereqs subRunData(targetDataLocation);

                subRunData.readData(inFile, "lowThresholdDet", j);

                // find the correct CSPrereqs to add this target's data to
                for(int k=0; (size_t)k<allData.size(); k++)
                {
                    if(allData[k].target.getName() == subRunData.target.getName())
                    {
                        // add subrun data to total
                        if(!allData[k].energyHisto)
                        {
                            // this is the first subrun to be added
                            allData[k].monitorCounts = subRunData.monitorCounts;
                            allData[k].energyHisto = (TH1I*)subRunData.energyHisto->Clone();
                            // prevent the cloned histogram from being closed
                            // when the subrun is closed
                            allData[k].energyHisto->SetDirectory(0);
                        }

                        else
                        {
                            allData[k] = allData[k] + subRunData;
                        }

                        break;
                    }

                    if((size_t)k+1==allData.size())
                    {
                        cerr << "Failed to find a CSPrereqs to add this subrun to." << endl;
                        continue;
                    }
                }
            }

            // Close the sub-run input files
            inFile->Close();
        }
    }

    string outFileName = dataLocation + "/total.root";

    vector<CrossSection> crossSections;

    cout << endl << "Total statistics over all runs: " << endl << endl;

    for(CSPrereqs p : allData)
    {
        long totalCounts = 0;
        for(int i=0; i<p.energyHisto->GetNbinsX(); i++)
        {
            totalCounts += p.energyHisto->GetBinContent(i);
        }

        cout << p.target.getName() << ": total events in energy histo = "
             << totalCounts << ", total monitor events = "
             << p.monitorCounts << endl;
        crossSections.push_back(calculateCS(outFileName, p, allData[0]));
    }

    /*string relativeFileName = dataLocation + "/relative.root";
    TFile* relativeFile = new TFile(relativeFileName.c_str(), "RECREATE");

    CrossSection relative = (crossSections[3]-crossSections[1])/(crossSections[3]+crossSections[1]);
    string relativeName = "#frac{#sigma_{" + allData[3].target.getName() +
                          "}-#sigma_{" + allData[1].target.getName() + 
                          "}}{#sigma_{" + allData[3].target.getName() +
                          "}+#sigma_{" + allData[1].target.getName() + "}}";

    relative.createCSGraph(relativeName.c_str());
    relativeFile->Close();
    */
}
