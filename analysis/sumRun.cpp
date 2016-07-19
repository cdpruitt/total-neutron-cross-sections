#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMath.h"

using namespace std;

const int NUMBER_OF_TARGETS = 6;
const int MAX_SUBRUN_NUMBER = 30;

const vector<string> targetNames = {"blankDPP", "shortCarbonDPP", "longCarbonDPP", "Sn112DPP", "NatSnDPP", "Sn124DPP"}; 
const vector<string> targetNamesWaveform = {"blankWaveform", "shortCarbonWaveform", "longCarbonWaveform", "Sn112Waveform", "NatSnWaveform", "Sn124Waveform"}; 

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

int main(int argc, char *argv[])
{
    string runNumber = argv[1];
    string driveName = argv[2];

    // create vectors for holding cross section data:
    // crossSections[target number]->at(subrun number)->at(data point)
    vector<vector<vector<double>*>*> crossSections;
    vector<vector<vector<double>*>*> crossSectionsWaveform;

    // create vectors for holding cross section average over all subruns:
    // crossSectionsAvg[target number]->at(data point)
    vector<vector<double>*> crossSectionsAvg;
    vector<vector<double>*> crossSectionsWaveformAvg;

    // create vector for holding the energies were the cross sections
    // were calculated
    vector<vector<vector<double>*>*> energies;
    vector<vector<vector<double>*>*> energiesWaveform;

    // prep vectors for filling
    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        crossSections.push_back(new vector<vector<double>*>);
        crossSectionsWaveform.push_back(new vector<vector<double>*>);

        crossSectionsAvg.push_back(new vector<double>);
        crossSectionsWaveformAvg.push_back(new vector<double>);

        energies.push_back(new vector<vector<double>*>);
        energiesWaveform.push_back(new vector<vector<double>*>);
    }

    TFile* infile;
    TFile* infileWaveform;

    // Loop over subruns to create average histos
    for(int i = 0; i<=MAX_SUBRUN_NUMBER; i++)
    {
        // Calculate the name of the next sub-run
        if(i < 10)
        {
            infile =  new TFile(Form("%s/analysis/run%s/run%s-000%i_cross-sections.root",driveName.c_str(),runNumber.c_str(),runNumber.c_str(),i));
            infileWaveform =  new TFile(Form("%s/analysis/run%s/run%s-000%i_waveform.root",driveName.c_str(),runNumber.c_str(),runNumber.c_str(),i));
        }

        else if(i < 100)
        {
            infile =  new TFile(Form("%s/analysis/run%s/run%s-00%i_cross-sections.root",driveName.c_str(),runNumber.c_str(),runNumber.c_str(),i));
            infileWaveform =  new TFile(Form("%s/analysis/run%s/run%s-00%i_waveform.root",driveName.c_str(),runNumber.c_str(),runNumber.c_str(),i));
        }

        else
        {
            cout << "Error: subrun number too large." << endl;
            infile->Close();
            infileWaveform->Close();
            exit(1);
        }

        // Attempt to open the sub-run
        if(!infile->IsOpen())
        {
            cout << "Can't open subrun " << runNumber << " " << i << endl;
            continue;
        }

        // Successfully found the sub-run of interest
        cout << "Adding run " << runNumber << " " << i <<endl;

        // Pull out the cross section data
        for(int j=0; j<targetNames.size(); j++)
        {
            TGraph * graph = (TGraph*)infile->Get(targetNames[j].c_str());
            energies[j]->push_back(new vector<double>);
            crossSections[j]->push_back(new vector<double>);
            extractGraphData(graph,energies[j]->back(),crossSections[j]->back());

            TGraph * graphWaveform = (TGraph*)infileWaveform->Get(targetNamesWaveform[j].c_str());
            energiesWaveform[j]->push_back(new vector<double>);
            crossSectionsWaveform[j]->push_back(new vector<double>);
            extractGraphData(graphWaveform,energiesWaveform[j]->back(),crossSectionsWaveform[j]->back());
        }

        infile->Close();
        infileWaveform->Close();
        // End of loop - move to next sub-run
    }

    // Create output file to contain summed histos
    stringstream outfileName;
    outfileName << driveName << "/analysis/run" << runNumber << "/" << "sum.root";
    TFile *outfile = new TFile(outfileName.str().c_str(),"RECREATE");

    // sum all subruns and compute average
    for(int i=1; i<energies.size(); i++)
    {
        crossSectionsAvg[i]->resize(energies[0]->at(0)->size());
        for(int j=0; j<energies[1]->size(); j++)
        {
            for(int k=0; k<energies[1]->at(0)->size(); k++)
            {
                crossSectionsAvg[i]->at(k) += crossSections[i]->at(j)->at(k);
            }
        }

        for(int k=0; k<energies[1]->at(0)->size(); k++)
        {
            crossSectionsAvg[i]->at(k) /= energies[1]->size();
        }

        // create new graphs to display the average
        TGraph* graph = new TGraph(energies[i]->at(0)->size(),&energies[i]->at(0)->at(0),&crossSectionsAvg[i]->at(0));
        graph->SetNameTitle(targetNames[i].c_str(),targetNames[i].c_str());
        graph->Write();
    }

    // sum all subruns and compute average
    for(int i=0; i<energiesWaveform.size(); i++)
    {
        crossSectionsWaveformAvg[i]->resize(energiesWaveform[0]->at(0)->size());
        for(int j=0; j<energiesWaveform[0]->size(); j++)
        {
            for(int k=0; k<energiesWaveform[0]->at(0)->size(); k++)
            {
                crossSectionsWaveformAvg[i]->at(k) += crossSectionsWaveform[i]->at(j)->at(k);
            }
        }

        for(int k=0; k<energiesWaveform[0]->at(0)->size(); k++)
        {
            crossSectionsWaveformAvg[i]->at(k) /= energiesWaveform[0]->size();
        }

        // create new graphs to display the average
        TGraph* graph = new TGraph(energiesWaveform[i]->at(0)->size(),&energiesWaveform[i]->at(0)->at(0),&crossSectionsWaveformAvg[i]->at(0));
        graph->SetNameTitle(targetNamesWaveform[i].c_str(),targetNamesWaveform[i].c_str());
        graph->Write();
    }

    outfile->Write();
    outfile->Close();
}
