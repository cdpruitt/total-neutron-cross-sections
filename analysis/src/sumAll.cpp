#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TLatex.h"

#include "../include/targetConstants.h"
#include "../include/dataSet.h"
#include "../include/dataPoint.h"

using namespace std;

const double ARTIFICIAL_OFFSET = 0.0;

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

TFile* infile;

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

void createRelativeCSPlot(DataSet set1, DataSet set2,
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
}

void createRelativeCSPlot(vector<double>* graph1Data, vector<double>* graph1Error,
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
}

DataSet scaleToLit(DataSet setToScale, DataSet expReference, DataSet litReference)
{
    DataSet scaledSet;
    int n = setToScale.getNumberOfPoints();
    for(int i=0; i<n; i++)
    {
        DataPoint p = litReference.getPoint(i)/expReference.getPoint(i);
        scaledSet.addPoint(p);
    }

    return scaledSet;
}

int main(int, char* argv[])
{
    // Find the total number of runs to read in
    string runNumber;
    ifstream runList;

    vector<vector<DataSet*>*> dataSets;

    // Loop through all listed runs and extract cross section data

    string targetName = argv[1];
    targetName = "../" + targetName + "/runsToSort.txt";

    runList.open(targetName);
    int i = 0;
    while (runList >> runNumber)
    {
        // Determine correct drive to find run
        string driveName;
        if (stoi(runNumber)<6)
        {
            driveName = "/data3/analysis/run";
        }
        else if (stoi(runNumber)>127 && stoi(runNumber)<160)
        {
            driveName = "/data2/analysis/run";
        }
        else if (stoi(runNumber)>159 && stoi(runNumber)<178)
        {
            driveName = "/data3/analysis/run";
        }
        else
        {
            cout << "Run directory outside bounds (runs 128-177) - check run number" << endl;
            exit(1);
        }

        // Open run

        TFile* infile;

        stringstream infileName;
        infileName << driveName << runNumber << "/" << "sum.root";
        infile = new TFile(infileName.str().c_str());

        if(!infile->IsOpen())
        {
            cout << "Can't open sum.root file of run" << runNumber << endl;
            exit(1);
        }

        cout << "Adding run" << runNumber << endl;

        dataSets.push_back(new vector<DataSet*>);

        // Pull out the cross section data
        for(int j=1; (size_t)j<targetNames.size(); j++)
        {
            TGraphErrors* graph = (TGraphErrors*)infile->Get(targetNames[j].c_str());
            
            dataSets.back()->push_back(new DataSet);
            extractGraphData(graph, *dataSets.back()->back());
        }

        // Close the sub-run input files
        infile->Close();

        // End of loop - move to next sub-run
        i++;
    }

    // read literature data for natural Sn
    //TFile *litData = new TFile("/data2/analysis/literatureData.root","READ");
    //litData->ls();
    //TGraphErrors *SnNatLitData = (TGraphErrors*)litData->Get("Natural Sn (n,tot)");

    // Create output file to contain averaged TGraphs of experimental cross section data
    stringstream outFileName;
    outFileName << "/data3/analysis/total.root";
    TFile *outFile = new TFile(outFileName.str().c_str(),"RECREATE");

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
    TFile *litDataFile = new TFile("/data2/analysis/literatureData.root","READ");
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

    createRelativeCSPlot(dataAverage[1],dataAverage[0],"long-to-short carbon relative CS",
                         "long/short carbon relative CS",
                         false);

    createRelativeCSPlot(dataAverage[0],CNatLitData,"short carbon/lit carbon relative CS",
                         "short carbon/lit carbon relative CS",
                         false);

    createRelativeCSPlot(dataAverage[0],CNatLitData,"long/short carbon relative CS",
                         "long/short carbon relative CS",
                         false);

    // create isotopic Sn plots, scaled using literature data for natural Sn
    DataSet Sn112Scaled = scaleToLit(dataAverage[2],dataAverage[3],SnNatLitData);
    DataSet Sn124Scaled = scaleToLit(dataAverage[4],dataAverage[3],SnNatLitData);

    Sn112Scaled.createPlot("Sn112Scaled");
    Sn124Scaled.createPlot("Sn124Scaled");

    /*for(int i=0; i<crossSectionsAvg[4].size(); i++)
    {
        DataPoint dataPoint(crossSectionsAvg[3][i]
        Sn112Scaled.addDataPoint(DataPoint(crossSectionsAvg[]))
    }*/

    // Write run histograms to sum.root
    outFile->Write();
    outFile->Close();
}
