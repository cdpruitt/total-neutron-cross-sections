#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"

using namespace std;

const int NUMBER_OF_TARGETS = 6;
const int MAX_SUBRUN_NUMBER = 30;
const double ARTIFICIAL_OFFSET = 0.0;

const vector<string> targetNames = {"blankDPP", "shortCarbonDPP", "longCarbonDPP", "Sn112DPP", "NatSnDPP", "Sn124DPP"}; 
const vector<string> targetNamesWaveform = {"blankWaveform", "shortCarbonWaveform", "longCarbonWaveform", "Sn112Waveform", "NatSnWaveform", "Sn124Waveform"}; 

// Set number of bins for relative cross-section histograms
const int NUM_RELATIVE_BINS = 50;

struct Plots
{
    vector<TGraph*> CSGraphs;
    vector<TGraph*> CSGraphsWaveform;

    vector<TGraph*> relativeCSGraphs;
    vector<TGraph*> relativeCSGraphsWaveform;
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

/*void getStatisticalError(TH1D* Sn112CSTotal, TH1D* Sn124CSTotal)
{
    // Get monCounts, raw target counts from _histos.root files
    // Calculate statistical error from histogram counts

    int limit = 25;

    ifstream runList("runsToSort.txt");

    string runDir;
    string outPath;

    while (runList >> runDir)
    {
        if (stoi(runDir)<6)
        {
            outPath = "/data3/analysis/run";
        }

        else if (stoi(runDir)>128 && stoi(runDir)<160)
        {
            outPath = "/data2/analysis/run";
        }

        else if (stoi(runDir)>159 && stoi(runDir)<178)
        {
            outPath = "/data3/analysis/run";
        }

        else
        {
            cout << "Run directory outside bounds (runs 128-177) - check run number" << endl;
            exit(1);
        }

        // keep track of which order the targets are in, based on which run number we're
        // sorting
        vector<int> order;

        if(stoi(runDir)<=5)
        {
            cout << "Neon run - stopping sort." << endl;
            exit(0);
        }

        if(stoi(runDir)<=151)
        {
            // blank, short carbon, long carbon, Sn112, NatSn, Sn124
            order.push_back(0);
            order.push_back(1);
            order.push_back(2);
            order.push_back(3);
            order.push_back(4);
            order.push_back(5);
        }

        else if(stoi(runDir)==152)
        {
            // blank, Sn112, NatSn, Sn124, short carbon, long carbon
            order.push_back(0);
            order.push_back(3);
            order.push_back(4);
            order.push_back(5);
            order.push_back(1);
            order.push_back(2);
        }

        else if(stoi(runDir)>=153 && stoi(runDir)<=168)
        {
            // blank, Sn112, NatSn, Sn124
            order.push_back(0);
            order.push_back(3);
            order.push_back(4);
            order.push_back(5);
        }

        else if(stoi(runDir)>=169 && stoi(runDir)<=180)
        {
            // blank, Sn112, NatSn, Sn124, short carbon
            order.push_back(0);
            order.push_back(3);
            order.push_back(4);
            order.push_back(5);
            order.push_back(1);
        }

        int pos112 = find(order.begin(), order.end(), 3) - order.begin();
        int pos124 = find(order.begin(), order.end(), 5) - order.begin();

        for(int segment = 0; segment<=limit; segment++)
        {
            // We need to form the proper name for the sub-run we want to open:
            if(segment < 10)
            {
                infile =  new TFile(Form("%s%s/run%s-000%i_histos.root",outPath.c_str(),runDir.c_str(),runDir.c_str(),segment));
            }

            else if(segment < 100)
            {
                infile =  new TFile(Form("%s%s/run%s-00%i_histos.root",outPath.c_str(),runDir.c_str(),runDir.c_str(),segment));
            }

            else
            {
                // There's some error in the sub-run file number; write outfile and exit
                cout << "Segment number too large!" << endl;
                exit(1);
            }

            // Attempt to open the sub-run to access its histograms
            if(!infile->IsOpen())
            {
                cout << "Can't open root file run " << runDir << " segment = " << segment << endl;

                break;
            }

            gDirectory->cd("/");
            gDirectory->GetDirectory("monitor")->cd();

            monCountsBlank += ((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(2);
            monCountsSn112 += ((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(2+pos112);
            monCountsSn124 += ((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(2+pos124);

            gDirectory->cd("/");
            gDirectory->GetDirectory("detS")->cd();

            TH1I* blank = ((TH1I*)gDirectory->Get("blankLog"));

            string name112 = "target" + to_string(pos112) + "Log";
            string name124 = "target" + to_string(pos124) + "Log";

            TH1I* Sn112 = ((TH1I*)gDirectory->Get(name112.c_str()));
            TH1I* Sn124 = ((TH1I*)gDirectory->Get(name124.c_str()));

            for(int i=0; i<blank->GetNbinsX(); i++)
            {
                detCountsBlank[i] += (blank->GetBinContent(i));
                detCountsSn112[i] += (Sn112->GetBinContent(i));
                detCountsSn124[i] += (Sn124->GetBinContent(i));
            }

            // Close the sub-run input files
            infile->Close();

            // End of loop - move to next sub-run
        }
    }

    *//*monCountsBlankError = pow(monCountsBlank,0.5);
    monCountsSn112Error = pow(monCountsSn112,0.5);
    monCountsSn124Error = pow(monCountsSn124,0.5);
    *//*

    cout << "monCountsBlank = " <<  monCountsBlank << endl;
    cout << "monCountsSn112 = " <<  monCountsSn112 << endl;
    cout << "monCountsSn124 = " <<  monCountsSn124 << endl;

    cout << "detCountsBlank[10] = " << detCountsBlank[10] << endl;
    cout << "detCountsSn112[10] = " << detCountsSn112[10] << endl;
    cout << "detCountsSn124[10] = " << detCountsSn124[10] << endl;


    for(int i=0; (size_t)i<totalSn112StatError.size(); i++)
    {
        totalSn112StatError[i] = pow(
                ((double)monCountsBlank/(pow(monCountsBlank,2))) +
                ((double)monCountsSn112/(pow(monCountsSn112,2))) +
                ((double)detCountsBlank[i]/(pow(detCountsBlank[i],2))) +
                ((double)detCountsSn112[i]/(pow(detCountsSn112[i],2)))
                ,0.5)*//*/exp(abs(Sn112CSTotal->GetBinContent(i)))*//*;

        totalSn124StatError[i] = pow(
                ((double)monCountsBlank/(pow(monCountsBlank,2))) +
                ((double)monCountsSn124/(pow(monCountsSn124,2))) +
                ((double)detCountsBlank[i]/(pow(detCountsBlank[i],2))) +
                ((double)detCountsSn124[i]/(pow(detCountsSn124[i],2)))
                ,0.5)*//*/exp(abs(Sn124CSTotal->GetBinContent(i)))*//*;

        *//*totalRelStatError[i] = pow(
                pow(totalSn112Error[i],2) +
                pow(totalSn124Error[i],2)
                ,0.5);
                *//*
    }

    cout << totalSn112StatError[10] << endl;

        *//*detCountsBlankpush_back(pow(detCountsBlank[i],0.5));
        detCountsSn112push_back(pow(detCountsSn112[i],0.5));
        detCountsSn124push_back(pow(detCountsSn124[i],0.5));
        *//*
}

void getSystemicError()
{
    // physical target data, listed in order:
    // {blank, Sn112, Natural Sn, Sn124, short carbon, long carbon} 

    // lengths of each target:
    double targetlength[6] = {0,1.37,2.74,1.365,1.370,1.370}; //cm
    double targetlengthUncertainty[6] = {0,0.01,0.01,0.005,0.005,0.005}; //cm

    // molar mass of each target:
    double targetMolMass[6] = {0,12.01,12.01,112,118.7,124}; //g/mol
    double targetMolMassUncertainty[6] = {0,0.01,0.01,0.01,0.01,0.01}; //g/mol

    // density of each target:
    double targetdensity[6] = {0,2.2,2.2,6.89,7.31,7.63}; //g/cm^3
    double targetdensityUncertainty[6] = {0,0.1,0.1,0.001,0.001,0.001}; //g/cm^3

    totalSn112SysError = pow(
            pow(targetlengthUncertainty[3]/targetlength[3],2) +
            pow(targetMolMassUncertainty[3]/targetMolMass[3],2) +
            pow(targetdensityUncertainty[3]/targetdensity[3],2)
            ,0.5);

    totalSn124SysError = pow(
            pow(targetlengthUncertainty[5]/targetlength[5],2) +
            pow(targetMolMassUncertainty[5]/targetMolMass[5],2) +
            pow(targetdensityUncertainty[5]/targetdensity[5],2)
            ,0.5);
}

void getTotalError()
{
    for(int i=0; (size_t)i<totalSn112Error.size(); i++)
    {
        totalSn112Error[i] = pow(
                pow(totalSn112StatError[i],2) +
                pow(totalSn112SysError,2)
                ,0.5);
        totalSn124Error[i] = pow(
                pow(totalSn124StatError[i],2) +
                pow(totalSn124SysError,2)
                ,0.5);
        totalRelError[i] = pow(
                pow(totalSn112Error[i],2) +
                pow(totalSn124Error[i],2)
                ,0.5);
    }

    cout << totalSn112Error[10] << endl;
    cout << totalSn124Error[10] << endl;
    cout << totalRelError[10] << endl;
}
*/


// Extract each point from a graph and store their positions in two vectors,
// xValues and yValues
void extractGraphData(
        TGraph* graph,
        vector<double>* xValues,
        vector<double>* yValues)
{
    int numPoints = graph->GetN();
    yValues->resize(numPoints);
    xValues->resize(numPoints);

    for(int k=0; k<numPoints; k++)
    {
        graph->GetPoint(k,xValues->at(k),yValues->at(k));
    }
}

int main()
{
    // Find the total number of runs to read in
    int totalRuns = 0;
    string runNumber;
    ifstream runList;

    runList.open("runsToSort.txt");
    while(runList >> runNumber)
    {
        totalRuns++;
    }
    cout << "Total runs in runlist: " << totalRuns << endl << endl;
    runList.close();

    // Create vectors for holding cross section data:
    // CrossSections[target number]->at(subrun number)->at(data point)
    vector<vector<vector<double>*>*> crossSections;
    vector<vector<vector<double>*>*> crossSectionsWaveform;

    // Create vectors for holding cross section average over all subruns:
    // CrossSectionsAvg[target number]->at(data point)
    vector<vector<double>*> crossSectionsAvg;
    vector<vector<double>*> crossSectionsWaveformAvg;

    // Create vector for holding the energies were the cross sections
    // were calculated
    vector<vector<vector<double>*>*> energies;
    vector<vector<vector<double>*>*> energiesWaveform;

    // Prep vectors for filling
    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        crossSections.push_back(new vector<vector<double>*>);
        crossSectionsWaveform.push_back(new vector<vector<double>*>);

        crossSectionsAvg.push_back(new vector<double>);
        crossSectionsWaveformAvg.push_back(new vector<double>);

        energies.push_back(new vector<vector<double>*>);
        energiesWaveform.push_back(new vector<vector<double>*>);
    }

    // Loop through all listed runs and extract cross section data
    runList.open("runsToSort.txt");
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

        // Pull out the cross section data
        for(int j=1; (size_t)j<targetNames.size(); j++)
        {
            TGraph * graph = (TGraph*)infile->Get(targetNames[j].c_str());
            energies[j]->push_back(new vector<double>);
            crossSections[j]->push_back(new vector<double>);
            extractGraphData(graph,energies[j]->back(),crossSections[j]->back());

            TGraph * graphWaveform = (TGraph*)infile->Get(targetNamesWaveform[j].c_str());
            energiesWaveform[j]->push_back(new vector<double>);
            crossSectionsWaveform[j]->push_back(new vector<double>);
            extractGraphData(graphWaveform,energiesWaveform[j]->back(),crossSectionsWaveform[j]->back());
        }

        // Close the sub-run input files
        infile->Close();

        // End of loop - move to next sub-run
        i++;
    }

    // Create output file to contain summed histos
    stringstream outfileName;
    outfileName << "/data3/analysis/total.root";
    TFile *outfile = new TFile(outfileName.str().c_str(),"RECREATE");

    // sum all subruns and compute average
    for(int i=1; (size_t)i<energies.size(); i++)
    {
        crossSectionsAvg[i]->resize(energies[1]->at(0)->size());
        for(int j=0; (size_t)j<energies[1]->size(); j++)
        {
            for(int k=0; (size_t)k<energies[1]->at(0)->size(); k++)
            {
                crossSectionsAvg[i]->at(k) += ARTIFICIAL_OFFSET+crossSections[i]->at(j)->at(k);
            }
        }

        for(int k=0; (size_t)k<energies[1]->at(0)->size(); k++)
        {
            crossSectionsAvg[i]->at(k) /= energies[1]->size();
        }

        // create new graphs to display the average
        TGraph* graph = new TGraph(energies[i]->at(0)->size(),&energies[i]->at(0)->at(0),&crossSectionsAvg[i]->at(0));
        graph->SetNameTitle(targetNames[i].c_str(),targetNames[i].c_str());
        graph->Write();
        plots.CSGraphs.push_back(graph);
    }

    // sum all subruns and compute average
    for(int i=1; (size_t)i<energiesWaveform.size(); i++)
    {
        crossSectionsWaveformAvg[i]->resize(energiesWaveform[1]->at(0)->size());
        for(int j=0; (size_t)j<energiesWaveform[1]->size(); j++)
        {
            for(int k=0; (size_t)k<energiesWaveform[1]->at(0)->size(); k++)
            {
                crossSectionsWaveformAvg[i]->at(k) += ARTIFICIAL_OFFSET+crossSectionsWaveform[i]->at(j)->at(k);
            }
        }

        for(int k=0; (size_t)k<energiesWaveform[1]->at(0)->size(); k++)
        {
            crossSectionsWaveformAvg[i]->at(k) /= energiesWaveform[1]->size();
        }

        // create new graphs to display the average
        TGraph* graph = new TGraph(energiesWaveform[i]->at(0)->size(),&energiesWaveform[i]->at(0)->at(0),&crossSectionsWaveformAvg[i]->at(0));
        graph->SetNameTitle(targetNamesWaveform[i].c_str(),targetNamesWaveform[i].c_str());
        graph->Write();
        plots.CSGraphs.push_back(graph);
    }

    // create relative cross-section plot for 112Sn/124Sn
    vector<double> Sn124PlusSn112CSData; 
    vector<double> Sn124MinusSn112CSData; 
    vector<double> Sn124Sn112relativeCSData; 

    Sn124PlusSn112CSData.resize(crossSectionsAvg[1]->size());
    Sn124MinusSn112CSData.resize(crossSectionsAvg[1]->size());
    Sn124Sn112relativeCSData.resize(crossSectionsAvg[1]->size());

    for(int i=1; (size_t)i<Sn124PlusSn112CSData.size(); i++)
    {
        Sn124PlusSn112CSData[i] = crossSectionsAvg[5]->at(i) + crossSectionsAvg[3]->at(i);
        Sn124MinusSn112CSData[i] = crossSectionsAvg[5]->at(i) - crossSectionsAvg[3]->at(i);
        Sn124Sn112relativeCSData[i] =
            (crossSectionsAvg[5]->at(i) - crossSectionsAvg[3]->at(i))/
            (crossSectionsAvg[5]->at(i) + crossSectionsAvg[3]->at(i));
    }

    TGraph *Sn124PlusSn112CS =
        new TGraph(Sn124PlusSn112CSData.size(),
                  &energies[1]->at(0)->at(0),
                  &Sn124PlusSn112CSData[0]);
    Sn124PlusSn112CS->Write();

    TGraph *Sn124MinusSn112CS =
        new TGraph(Sn124MinusSn112CSData.size(),
                  &energies[1]->at(0)->at(0),
                  &Sn124MinusSn112CSData[0]);
    Sn124MinusSn112CS->Write();

    TGraph *Sn124Sn112relativeCS =
        new TGraph(Sn124Sn112relativeCSData.size(),
                &energies[1]->at(0)->at(0),
                &Sn124Sn112relativeCSData[0]);
    Sn124Sn112relativeCS->SetTitle("#frac{#sigma_{^{124}Sn}-#sigma_{^{112}Sn}}{#sigma_{^{124}Sn}+#sigma_{^{112}Sn}}");
    Sn124Sn112relativeCS->Write();

    //getStatisticalError(Sn112CSTotalLog, Sn124CSTotalLog);
    //getSystemicError();
    //getTotalError();

    /*vector<double> energy;
    vector<double> xsection;
    vector<double> error;
    */

    // Write run histograms to sum.root
    outfile->Write();
    outfile->Close();
}
