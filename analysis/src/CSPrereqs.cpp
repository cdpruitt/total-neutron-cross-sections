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

void CSPrereqs::getHisto(TFile* histoFile, string directory, string name)
{
    // for histos
    string energyHistoName = name + "Energy";

    histoFile->cd(directory.c_str());
    energyHisto = ((TH1D*)gDirectory->Get(energyHistoName.c_str()));

    if(!energyHisto)
    {
        cerr << "Error: Failed to open histogram \"" << energyHistoName
            << "\" in file " << histoFile
            << " (in getHisto)." << endl;
        return;
    }

    string TOFHistoName = name + "TOFCorrected";
    TOFHisto = ((TH1D*)gDirectory->Get(TOFHistoName.c_str()));

    if(!TOFHisto)
    {
        cerr << "Error: Failed to open histogram \"" << TOFHistoName
            << "\" in file " << histoFile
            << " (in getHisto)." << endl;
        return;
    }
}

// for histos
void CSPrereqs::getMonitorCounts(TFile* histoFile, string directory, int targetPosition)
{
    TDirectory* dir = (TDirectory*)histoFile->Get(directory.c_str());
    if(!dir)
    {
        cerr << "Error: Failed to open directory " << directory << " in file " << histoFile->GetName()
            << " (in getMonitorCounts)." << endl;
        return;
    }

    TH1D* targetPosHisto = (TH1D*)dir->Get("targetPosH");
    if(!targetPosHisto)
    {
        cerr << "Error: Failed to open histogram \"targetPosH\" in file " << histoFile
            << " (in getMonitorCounts)." << endl;
        return;
    }

    monitorCounts = targetPosHisto->GetBinContent(targetPosition+1);

    if(monitorCounts<=0)
    {
        cerr << "Error: for target position " << config.target.TARGET_ORDER[targetPosition] << ", monitor counts read was <=0 (in getMonitorCounts)." << endl;
        return;
    }
}

void CSPrereqs::getMacroNumber(TFile* histoFile, TFile* monitorFile, string directory, string goodMacroHistoName, string macroHistoName)
{
    TDirectory* dir = (TDirectory*)histoFile->Get(directory.c_str());
    if(!dir)
    {
        cerr << "Error: Failed to open directory " << directory << " in file " << histoFile->GetName()
            << " (in getMonitorCounts)." << endl;
        return;
    }

    TDirectory* monitorDir = (TDirectory*)monitorFile->Get(directory.c_str());
    if(!monitorDir)
    {
        cerr << "Error: Failed to open directory " << directory << " in file " << monitorFile->GetName()
            << " (in getMonitorCounts)." << endl;
        return;
    }

    TH1D* goodMacrosH = (TH1D*)dir->Get(goodMacroHistoName.c_str());
    if(!goodMacrosH)
    {
        cerr << "Error: Failed to open histogram \"" << goodMacroHistoName
            << "\" in file " << histoFile->GetName() << " (in getMacroRatio)." << endl;
        return;
    }

    unsigned int numberOfBins = goodMacrosH->GetNbinsX();

    unsigned int goodMacroCount = 0;

    for(int j=1; j<=numberOfBins; j++)
    {
        if(goodMacrosH->GetBinContent(j)>0)
        {
            goodMacroNumber++;
        }
    }

    TH1D* macroHisto = (TH1D*)monitorDir->Get(macroHistoName.c_str());
    if(!macroHisto)
    {
        cerr << "Error: Failed to open histogram \"" << macroHistoName
            << "\" in file " << histoFile->GetName() << " (in getMacroRatio)." << endl;
        return;
    }

    numberOfBins = macroHisto->GetNbinsX();

    unsigned int macroCount = 0;

    for(int j=1; j<numberOfBins; j++)
    {
        if(macroHisto->GetBinContent(j)>0)
        {
            totalMacroNumber++;
        }
    }

    if(goodMacroNumber<=0 || totalMacroNumber<=0)
    {
        cerr << "Error: total and/or good macro number cannot be <= 0." << endl;
    }
}

// readData for histos
void CSPrereqs::readEnergyData(TFile* histoFile, string directory, int targetPosition)
{
    // Find deadtime-corrected energy histo for this target
    string histoName = config.target.TARGET_ORDER[targetPosition];
    getHisto(histoFile, directory, histoName);
}

void CSPrereqs::readMonitorData(TFile* histoFile, string monitorDirectory, int targetPosition)
{
    // Find monitor histo for this target
    getMonitorCounts(histoFile, monitorDirectory, targetPosition);
}

void CSPrereqs::readMacroData(TFile* macroFile, TFile* monitorFile, string detectorName, int targetPosition)
{
    string goodMacroHistoName = config.target.TARGET_ORDER[targetPosition] + "GoodMacros";
    string macroHistoName = config.target.TARGET_ORDER[targetPosition] + "MacroNumber";

    getMacroNumber(macroFile, monitorFile, detectorName, goodMacroHistoName, macroHistoName);
}

CSPrereqs operator+(CSPrereqs& augend, CSPrereqs& addend)
{
    if(augend.target.getName() != addend.target.getName())
    {
        cerr << "Error: tried to added CSPrereqs of different targets." << endl;
        exit(1);
    }

    augend.energyHisto->Add(addend.energyHisto);
    augend.TOFHisto->Add(addend.TOFHisto);
    augend.monitorCounts += addend.monitorCounts;
    augend.goodMacroNumber += addend.goodMacroNumber;
    augend.totalMacroNumber += addend.totalMacroNumber;

    return augend;
}

CSPrereqs::CSPrereqs(Target t)
{
    target = t;
    TOFHisto = new TH1D("","",config.plot.TOF_BINS,config.plot.TOF_LOWER_BOUND,config.plot.TOF_UPPER_BOUND),target.getName();
    TOFHisto->SetDirectory(0);
    energyHisto = timeBinsToRKEBins(TOFHisto,target.getName());
    energyHisto->SetDirectory(0);

    monitorCounts = 0;
    goodMacroNumber = 0;
    totalMacroNumber = 0;
}

CSPrereqs::CSPrereqs(string targetDataLocation) : target(targetDataLocation)
{
    TOFHisto = new TH1D("","",config.plot.TOF_BINS,config.plot.TOF_LOWER_BOUND,config.plot.TOF_RANGE),target.getName();
    TOFHisto->SetDirectory(0);
    energyHisto = timeBinsToRKEBins(TOFHisto,target.getName());
    energyHisto->SetDirectory(0);

    monitorCounts = 0;
    goodMacroNumber = 0;
    totalMacroNumber = 0;
}

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

int readSubRun(vector<CSPrereqs>& allCSPrereqs, string expName, int runNumber, int subRun, string detectorName, string dataLocation)
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
        return 0;
    }

    // open subrun
    stringstream monitorFileName;
    monitorFileName << dataLocation << "/" << runNumber << "/"
        << subRunFormatted.str() << "/histos.root";
    ifstream f(monitorFileName.str());
    if(!f.good())
    {
        // failed to open this sub-run - skip to the next one
        cerr << "Couldn't open " << monitorFileName.str() << "; continuing.\r";
        fflush(stdout);
        return 0;
    }

    TFile* monitorFile = new TFile(monitorFileName.str().c_str(),"READ");

    stringstream energyFileName;
    energyFileName << dataLocation << "/" << runNumber << "/"
        << subRunFormatted.str() << "/energy.root";
    ifstream g(energyFileName.str());
    if(!g.good())
    {
        // failed to open this sub-run - skip to the next one
        cerr << "Couldn't open " << energyFileName.str() << "; continuing.\r";
        fflush(stdout);
        return 0;
    }

    g.close();

    TFile* energyFile = new TFile(energyFileName.str().c_str(),"READ");

    stringstream macroFileName;
    macroFileName << dataLocation << "/" << runNumber << "/"
        << subRunFormatted.str() << "/gatedHistos.root";
    ifstream h(macroFileName.str());
    if(!h.good())
    {
        // failed to open this sub-run - skip to the next one
        cerr << "Couldn't open " << macroFileName.str() << "; continuing.\r";
        fflush(stdout);
        return 0;
    }

    h.close();

    TFile* macroFile = new TFile(macroFileName.str().c_str(),"READ");

    // get target order for this run
    vector<string> targetOrder = getTargetOrder(expName, runNumber);

    cout << "Adding " << runNumber << ", subrun " << subRun << endl;

    // Loop through all target positions in this subrun
    for(int j=0; (size_t)j<targetOrder.size(); j++)
    {
        // pull data needed for CS calculation from subrun 
        string targetDataLocation = "../" + expName + "/targetData/" + targetOrder[j] + ".txt";
        CSPrereqs subRunData(targetDataLocation);

        subRunData.readEnergyData(energyFile, detectorName, j);
        subRunData.readMonitorData(monitorFile, "monitor", j);
        subRunData.readMacroData(macroFile, monitorFile, detectorName, j);

        // find the correct CSPrereqs to add this target's data to
        for(CSPrereqs& csp : allCSPrereqs)
        {
            if(csp.target.getName() == subRunData.target.getName())
            {
                // add subrun data to total
                csp = csp + subRunData;
            }
        }
    }

    // Close the sub-run input files
    energyFile->Close();
    monitorFile->Close();
    macroFile->Close();

    return 0;
}
