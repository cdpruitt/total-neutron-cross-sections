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
#include "../include/targetConstants.h"
#include "../include/dataSet.h"
#include "../include/dataPoint.h"
#include "../include/CSPrereqs.h"
#include "../include/analysisConstants.h"
#include "../include/crossSection.h"

using namespace std;

const int MAX_SUBRUN_NUMBER = 50;

struct Plots
{
    vector<TGraphErrors*> CSGraphs;
    vector<TGraphErrors*> CSGraphsWaveform;

    vector<TGraphErrors*> relativeCSGraphs;
    vector<TGraphErrors*> relativeCSGraphsWaveform;
} plots;

long monCountsBlank;
long monCountsSn124;
long monCountsSn112;

double monCountsBlankError;
double monCountsSn124Error;
double monCountsSn112Error;

vector<long> detCountsBlank (200);
vector<long> detCountsSn124 (200);
vector<long> detCountsSn112 (200);

vector<double> detCountsBlankError;
vector<double> detCountsSn124Error;
vector<double> detCountsSn112Error;

vector<double> totalSn112StatError (200);
vector<double> totalSn124StatError (200);
vector<double> totalRelStatError (200);

double totalSn112SysError;
double totalSn124SysError;
double totalRelSysError;

vector<double> totalSn112Error (200);
vector<double> totalSn124Error (200);
vector<double> totalRelError (200);

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

/*void createRelativeCSPlot(DataSet set1, DataSet set2,
                          string name, string title, bool listRMS)
{
    // the relative difference plots are defined as:
    // (graph1Data-graph2Data)/(graph1Data+graph2Data)

    // define vectors to hold the data points and error for the numerator
    // and denominator of the relative difference

    DataSet summed = set1+set2;
    DataSet difference = set1-set2;
    DataSet relative = difference/summed;

    TGraphErrors *summedGraph =
        new TGraphErrors(summed.getNumberOfPoints(),
                  &summed.getXValues()[0],
                  &summed.getYValues()[0],
                  &summed.getXErrors()[0],
                  &summed.getYErrors()[0]);
    summedGraph->SetName("sumCS");
    summedGraph->Write();

    TGraphErrors *differenceGraph =
        new TGraphErrors(difference.getNumberOfPoints(),
                  &difference.getXValues()[0],
                  &difference.getYValues()[0],
                  &difference.getXErrors()[0],
                  &difference.getYErrors()[0]);
    differenceGraph->SetName("diffCS");
    differenceGraph->Write();

    TGraphErrors *relativeGraph =
        new TGraphErrors(relative.getNumberOfPoints(),
                  &relative.getXValues()[0],
                  &relative.getYValues()[0],
                  &relative.getXErrors()[0],
                  &relative.getYErrors()[0]);
    relativeGraph->SetName("relativeCS");
    relativeGraph->Write();

    relativeGraph->SetLineColor(kBlue);
    relativeGraph->SetLineWidth(2);
    relativeGraph->Draw();

    // calculate RMS value
    if(listRMS)
    {
        double rms = calculateRMS(set1.getYValues(),set2.getYValues());

        // add RMS value to plot
        stringstream rmsText;
        rmsText << "rms = " << rms;

        TLatex* latex = new TLatex(0.3,0.2,rmsText.str().c_str());
        latex->SetNDC();
        latex->SetTextSize(0.05);
        relativeGraph->GetListOfFunctions()->Add(latex);
    }
    
    relativeGraph->Write();
}*/

/*void createRelativeCSPlot(vector<double>* graph1Data, vector<double>* graph1Error,
                          vector<double>* graph2Data, vector<double>* graph2Error,
                          vector<double>* energyData ,vector<double>* energyError,
                          string name, string title, bool listRMS)
{
    // the relative difference plots are defined as:
    // (graph1Data-graph2Data)/(graph1Data+graph2Data)

    // define vectors to hold the data points and error for the numerator
    // and denominator of the relative difference
    vector<double> summedGraphData; 
    vector<double> summedGraphError; 

    vector<double> differenceGraphData; 
    vector<double> differenceGraphError; 

    vector<double> relativeGraphData; 
    vector<double> relativeGraphError; 

    summedGraphData.resize(graph1Data->size());
    summedGraphError.resize(graph1Data->size());

    differenceGraphData.resize(graph1Data->size());
    differenceGraphError.resize(graph1Data->size());

    relativeGraphError.resize(graph1Data->size());
    relativeGraphData.resize(graph1Data->size());

    for(int i=1; (size_t)i<summedGraphData.size(); i++)
    {
        summedGraphData[i] = graph1Data->at(i) + graph2Data->at(i);
        differenceGraphData[i] = graph1Data->at(i) - graph2Data->at(i);
        relativeGraphData[i] =
            (graph1Data->at(i) - graph2Data->at(i))/
            (graph1Data->at(i) + graph2Data->at(i));
    }

    for(int i=1; (size_t)i<summedGraphError.size(); i++)
    {
        summedGraphError[i] = pow(pow(graph1Error->at(i),2) + pow(graph2Error->at(i),2),0.5);
        differenceGraphError[i] = summedGraphError[i];
        relativeGraphError[i] = relativeGraphData[i]*
                                       pow(
                                          pow(summedGraphError[i]/summedGraphData[i],2) +
                                          pow(differenceGraphError[i]/differenceGraphData[i],2),
                                          0.5);
        //graph1Datagraph2DatarelativeCSError[i] = differenceGraphError[i]*summedGraphData[i]+summedGraphError[i]*differenceGraphData[i];
    }

    TGraphErrors *summedGraph =
        new TGraphErrors(summedGraphData.size(),
                  &energyData->at(0),
                  &summedGraphData[0],
                  &energyError->at(0),
                  &summedGraphError[0]);
    summedGraph->SetName("sumCS");
    summedGraph->Write();

    TGraphErrors *differenceGraph =
        new TGraphErrors(differenceGraphData.size(),
                  &energyData->at(0),
                  &differenceGraphData[0],
                  &energyError->at(0),
                  &differenceGraphError[0]);
    differenceGraph->SetName("differenceCS");
    differenceGraph->Write();

    TGraphErrors *relativeGraph =
        new TGraphErrors(relativeGraphData.size(),
                &energyData->at(0),
                &relativeGraphData[0],
                &energyError->at(0),
                &relativeGraphError[0]);
    relativeGraph->SetTitle(title.c_str());
    relativeGraph->SetName(name.c_str());

    relativeGraph->SetLineColor(kBlue);
    relativeGraph->SetLineWidth(2);
    relativeGraph->Draw();

    // calculate RMS value
    if(listRMS)
    {
        double rms = calculateRMS(*graph1Data,*graph2Data);

        // add RMS value to plot
        stringstream rmsText;
        rmsText << "rms = " << rms;

        TLatex* latex = new TLatex(0.3,0.2,rmsText.str().c_str());
        latex->SetNDC();
        latex->SetTextSize(0.05);
        relativeGraph->GetListOfFunctions()->Add(latex);
    }
    
    //relativeGraph->Write();
}*/

/*DataSet scaleToLit(DataSet setToScale, DataSet expReference, DataSet litReference)
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
}*/

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

            cout << "Adding " << runNumber << ", subrun " << subRun << endl;

            // Loop through all target positions in this subrun
            for(int j=0; (size_t)j<targetOrder.size(); j++)
            {
                // pull data needed for CS calculation from subrun 
                string targetDataLocation = "../" + expName + "/targetData/" + targetOrder[j] + ".txt";
                CSPrereqs subRunData(targetDataLocation);

                subRunData.readData(inFile, get<1>(channelMap[5]), j);

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

    cout << "Total statistics over all runs: " << endl << endl;

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
        calculateCS(outFileName, p, allData[0]);
    }

    // read literature data for natural Sn
    //TFile *litData = new TFile("/data2/analysis/literatureData.root","READ");
    //litData->ls();
    //TGraphErrors *SnNatLitData = (TGraphErrors*)litData->Get("Natural Sn (n,tot)");

    /*

    vector<DataSet> dataAverage;
    for(int i=0; i<dataSets[0]->size(); i++)
    {
        dataAverage.push_back(*dataSets[0]->at(i));
    }

    // sum all subruns and compute average
    for(int i=0; (size_t)i<dataSets[0]->size(); i++)
    {
        for(int j=1; (size_t)j<dataSets.size(); j++)
        {
            DataSet d = *dataSets[j]->at(i);
            dataAverage[i] = dataAverage[i].merge(d);
        }

        TGraphErrors* graph = dataAverage[i].createPlot(targetNames[i+1]);
        plots.CSGraphs.push_back(graph);
    }
    */

    /*

    //getStatisticalError(Sn112CSTotalLog, Sn124CSTotalLog);
    //getSystemicError();
    //getTotalError();

    vector<double> energy;
    vector<double> xsection;
    vector<double> error;
    */

    // produce relative cross section with isotopic cross sections
    // adjusted by matching to the literature SnNat data

    // read literature data
/*    TFile *litDataFile = new TFile("/data2/analysis/literatureData.root","READ");
    TGraphErrors *SnNatLitGraph = (TGraphErrors*)litDataFile->Get("Natural Sn (n,tot)");
    TGraphErrors *CNatLitGraph = (TGraphErrors*)litDataFile->Get("Natural carbon (n,tot)");

    // extract literature data on natural Sn and natural carbon
    DataSet SnNatLitData;
    DataSet CNatLitData;

    vector<double> energies;
    for(int i=0; (size_t)i<dataSets[0]->at(0)->getNumberOfPoints(); i++)
    {
        energies.push_back(dataSets[0]->at(0)->getPoint(i).getXValue());
    }

    for(int i=0; (size_t)i<energies.size(); i++)
    {
        DataPoint d(energies[i],
                    0,
                    SnNatLitGraph->Eval(energies[i]),
                    SnNatLitGraph->GetErrorY(energies[i]));
        SnNatLitData.addPoint(d);

        DataPoint d2(energies[i],
                     0,
                     CNatLitGraph->Eval(energies[i]),
                     CNatLitGraph->GetErrorY(energies[i]));

        CNatLitData.addPoint(d2);
    }

    litDataFile->Close();
    outFile->cd();

    createRelativeCSPlot(dataAverage[4],dataAverage[2],"Sn relative CS",
                         "#frac{#sigma_{^{124}Sn}-#sigma_{^{112}Sn}}{#sigma_{^{124}Sn}+#sigma_{^{112}Sn}}",
                         false);
                         */

    /*createRelativeCSPlot(dataAverage[1],dataAverage[0],"long-to-short carbon relative CS",
                         "long/short carbon relative CS",
                         false);

    createRelativeCSPlot(dataAverage[0],CNatLitData,"short carbon/lit carbon relative CS",
                         "short carbon/lit carbon relative CS",
                         false);

    createRelativeCSPlot(dataAverage[0],CNatLitData,"long/short carbon relative CS",
                         "long/short carbon relative CS",
                         false);
    */
/*
    // create isotopic Sn plots, scaled using literature data for natural Sn
    DataSet Sn112Scaled = scaleToLit(dataAverage[2],dataAverage[3],SnNatLitData);
    DataSet Sn124Scaled = scaleToLit(dataAverage[4],dataAverage[3],SnNatLitData);

    Sn112Scaled.createPlot("Sn112Scaled");
    Sn124Scaled.createPlot("Sn124Scaled");
    */

    /*for(int i=0; i<crossSectionsAvg[4].size(); i++)
    {
        DataPoint dataPoint(crossSectionsAvg[3][i]
        Sn112Scaled.addDataPoint(DataPoint(crossSectionsAvg[]))
    }*/
}
