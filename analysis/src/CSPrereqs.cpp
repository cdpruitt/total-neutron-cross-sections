#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include "TFile.h"
#include "TH1D.h"
#include "TDirectory.h"

#include "../include/target.h"
#include "../include/plots.h"
#include "../include/CSPrereqs.h"
#include "../include/config.h"
#include "../include/experiment.h"

using namespace std;

extern Config config;

int CSPrereqs::readTOFHisto(TFile* histoFile, string directory, string targetName)
{
    // for histos
    histoFile->cd(directory.c_str());

    string TOFHistoName = targetName + "TOF";
    TH1D* TOFHistoTemp = (TH1D*)gDirectory->Get(TOFHistoName.c_str());

    if(!TOFHistoTemp)
    {
        cerr << "Error: Failed to open histogram \"" << TOFHistoName
            << "\" in file " << histoFile
            << " (in getHisto)." << endl;
        return 1;
    }

    TOFHisto = (TH1D*)TOFHistoTemp->Clone();
    TOFHisto->SetDirectory(0);

    return 0;
}

// for histos
int CSPrereqs::readMonitorCounts(TFile* histoFile, string directory, string targetName)
{
    TDirectory* dir = (TDirectory*)histoFile->Get(directory.c_str());
    if(!dir)
    {
        cerr << "Error: Failed to open directory " << directory << " in file " << histoFile->GetName()
            << " (in readMonitorCounts)." << endl;
        return 1;
    }

    string histoName = targetName + "GoodMacros";
    TH1D* histo = (TH1D*)dir->Get(histoName.c_str());
    if(!histo)
    {
        cerr << "Error: Failed to open histogram " << histoName << " in file " << histoFile->GetName()
            << " (in readMonitorCounts)." << endl;
        return 1;
    }

    int numberOfBins = histo->GetNbinsX();

    double totalMonitorCounts = 0;
    int totalGoodMonitorMacros = 0;

    for(int i=1; i<=numberOfBins; i++)
    {
        double binContent = histo->GetBinContent(i);

        if(binContent<=0.5)
        {
            continue;
        }

        totalMonitorCounts += binContent;
        totalGoodMonitorMacros++;
    }

    monitorCounts = totalMonitorCounts;
    goodMacroNumber = totalGoodMonitorMacros;

    return 0;
}

int CSPrereqs::readMacroData(TFile* histoFile, string directory, string totalMacrosHistoName)
{
    TDirectory* dir = (TDirectory*)histoFile->Get(directory.c_str());
    if(!dir)
    {
        cerr << "Error: Failed to open directory " << directory << " in file " << histoFile->GetName()
            << " (in getMacroNumber)." << endl;
        return 1;
    }

    /*TDirectory* monitorDir = (TDirectory*)monitorFile->Get(directory.c_str());
    if(!monitorDir)
    {
        cerr << "Error: Failed to open directory " << directory << " in file " << monitorFile->GetName()
            << " (in getMacroNumber)." << endl;
        return;
    }*/

    string name = totalMacrosHistoName + "MacroNumber";
    TH1D* totalMacrosHisto = (TH1D*)dir->Get(name.c_str());
    if(!totalMacrosHisto)
    {
        cerr << "Error: Failed to open histogram \"" << totalMacrosHistoName
            << "\" in file " << histoFile->GetName() << " (in getMacroRatio)." << endl;
        return 1;
    }

    int numberOfBins = totalMacrosHisto->GetNbinsX();

    for(int j=1; j<=numberOfBins; j++)
    {
        if(totalMacrosHisto->GetBinContent(j)>0)
        {
            totalMacroNumber++;
        }
    }

    return 0;
}

int CSPrereqs::readEventData(TFile* histoFile, string directoryName, string targetName)
{
    string eventDataHistoName = targetName + "TOF";

    TDirectory* dir = (TDirectory*)histoFile->Get(directoryName.c_str());
    if(!dir)
    {
        cerr << "Error: Failed to open directory " << directoryName << " in file " << histoFile->GetName()
            << " (in readEventData)." << endl;
        return 1;
    }

    TH1D* eventDataHisto = (TH1D*)dir->Get(eventDataHistoName.c_str());
    if(!eventDataHisto)
    {
        cerr << "Error: Failed to open histogram \"" << eventDataHistoName
            << "\" in file " << histoFile->GetName() << " (in readEventData)." << endl;
        return 1;
    }

   totalEventNumber = eventDataHisto->GetEntries();

   return 0;
}

int CSPrereqs::readUncorrectedTOFHisto(TFile* histoFile, string directoryName, string targetName)
{
    string eventDataHistoName = targetName + "TOF";

    TDirectory* dir = (TDirectory*)histoFile->Get(directoryName.c_str());
    if(!dir)
    {
        cerr << "Error: Failed to open directory " << directoryName << " in file " << histoFile->GetName()
            << " (in readEventData)." << endl;
        return 1;
    }

    TH1D* eventDataHistoTemp = (TH1D*)dir->Get(eventDataHistoName.c_str());
    if(!eventDataHistoTemp)
    {
        cerr << "Error: Failed to open histogram \"" << eventDataHistoName
            << "\" in file " << histoFile->GetName() << " (in readEventData)." << endl;
        return 1;
    }

    string clonedHistoName = eventDataHistoName + "uncorrected";
    uncorrectedTOFHisto = (TH1D*)eventDataHistoTemp->Clone(clonedHistoName.c_str());
    uncorrectedTOFHisto->SetDirectory(0);

   return 0;
}

CSPrereqs operator+(CSPrereqs& augend, CSPrereqs& addend)
{
    if(augend.target.getName() != addend.target.getName())
    {
        cerr << "Error: tried to added CSPrereqs of different targets." << endl;
        exit(1);
    }

    if(!addend.TOFHisto || !addend.monitorCounts || !addend.goodMacroNumber || !addend.uncorrectedTOFHisto)
    {
        return augend;
    }

    augend.TOFHisto->Add(addend.TOFHisto);
    augend.uncorrectedTOFHisto->Add(addend.uncorrectedTOFHisto);

    augend.monitorCounts += addend.monitorCounts;
    augend.goodMacroNumber += addend.goodMacroNumber;
    augend.totalMacroNumber += addend.totalMacroNumber;
    augend.totalEventNumber += addend.totalEventNumber;

    return augend;
}

CSPrereqs::CSPrereqs(Target t)
{
    target = t;
    TOFHisto = new TH1D("","",config.plot.TOF_BINS,config.plot.TOF_LOWER_BOUND,config.plot.TOF_UPPER_BOUND);
    TOFHisto->SetDirectory(0);

    uncorrectedTOFHisto = new TH1D("","",config.plot.TOF_BINS,config.plot.TOF_LOWER_BOUND,config.plot.TOF_UPPER_BOUND);
    uncorrectedTOFHisto->SetDirectory(0);

    monitorCounts = 0;
    goodMacroNumber = 0;
    totalMacroNumber = 0;
}

CSPrereqs::CSPrereqs(string targetDataLocation) : target(targetDataLocation)
{
    TOFHisto = new TH1D("","",config.plot.TOF_BINS,config.plot.TOF_LOWER_BOUND,config.plot.TOF_RANGE);
    TOFHisto->SetDirectory(0);

    uncorrectedTOFHisto = new TH1D("","",config.plot.TOF_BINS,config.plot.TOF_LOWER_BOUND,config.plot.TOF_UPPER_BOUND);
    uncorrectedTOFHisto->SetDirectory(0);

    monitorCounts = 0;
    goodMacroNumber = 0;
    totalMacroNumber = 0;
}

// Extract each point from a graph and store their positions in two vectors,
// xValues and yValues
void extractGraphData(
        TGraphAsymmErrors* graph,
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
        TGraphAsymmErrors* graph,
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

int readTargetData(vector<CSPrereqs>& allCSPrereqs, string expName)
{
    // check to see if data with this target has already been read in
    // (from previous runs). If not, create a new CSPrereq to hold the new
    // target's data
    for(string targetName : config.target.TARGET_ORDER)
    {
        bool CSPAlreadyExists = false;

        for(CSPrereqs csp : allCSPrereqs)
        {
            if(csp.target.getName() == targetName)
            {
                // CSPrereq for this target already exists
                CSPAlreadyExists = true;
                break;
            }
        }

        if(CSPAlreadyExists)
        {
            continue;
        }

        string targetDataLocation = "../" + expName + "/targetData/" + targetName + ".txt";
        allCSPrereqs.push_back(CSPrereqs(targetDataLocation));
    }

    return 0;
}

int readSubRun(CSPrereqs& subRunData, string expName, int runNumber, int subRun, string detectorName, string dataLocation)
{
    stringstream subRunFormatted;
    subRunFormatted << setfill('0') << setw(4) << subRun;

    // Skip subruns on the blacklist
    string blacklistDataLocation = "../" + expName + "/blacklist.txt";
    ifstream blacklist(blacklistDataLocation);

    if(!blacklist.good())
    {
        cerr << "Error: couldn't find blacklist in " << blacklistDataLocation << endl;
        return 1;
    }

    bool skip = false;
    string str;

    while(getline(blacklist,str))
    {
        string runSubrun = to_string(runNumber) + "-" + subRunFormatted.str();
        if(runSubrun == str)
        {
            cout << "Found sub-run " << runSubrun << " on blacklist; skipping..." << endl;

            skip = true;
            break;
        }
    }

    if(skip)
    {
        return 1;
    }

    // open file containing deadtime-corrected TOF histograms
    stringstream TOFFileName;
    TOFFileName << dataLocation << "/" << runNumber << "/"
        << subRunFormatted.str() << "/correctedHistos.root";
    ifstream g(TOFFileName.str());
    if(!g.good())
    {
        // failed to open this sub-run - skip to the next one
        cerr << "Couldn't open " << TOFFileName.str() << "; continuing.\r";
        fflush(stdout);
        return 1;
    }

    g.close();

    TFile* TOFFile = new TFile(TOFFileName.str().c_str(),"READ");

    // open file containing macropulse number data
    stringstream macroFileName;
    macroFileName << dataLocation << "/" << runNumber << "/"
        << subRunFormatted.str() << "/gatedHistos.root";
    ifstream h(macroFileName.str());
    if(!h.good())
    {
        // failed to open this sub-run - skip to the next one
        cerr << "Couldn't open " << macroFileName.str() << "; continuing.\r";
        fflush(stdout);
        return 1;
    }

    h.close();

    TFile* macroFile = new TFile(macroFileName.str().c_str(),"READ");

    // pull data needed for CS calculation from subrun 
    string targetName = subRunData.target.getName();

    if(subRunData.readTOFHisto(TOFFile, detectorName, targetName)
            || subRunData.readMonitorCounts(macroFile, "monitor", targetName)
            || subRunData.readEventData(macroFile, detectorName, targetName)
            || subRunData.readUncorrectedTOFHisto(macroFile, detectorName, targetName))
    {
        // error: failed to read one of the essential cross section quantities
        // for this target.

        // Close the sub-run input file
        TOFFile->Close();
        macroFile->Close();

        return 1;
    }

    // successfully read all subrun quantities required for CS calculation
    TOFFile->Close();
    macroFile->Close();

    return 0;
}
