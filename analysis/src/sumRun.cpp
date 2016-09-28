#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TMath.h"

using namespace std;

const int NUMBER_OF_TARGETS = 6;
const int MAX_SUBRUN_NUMBER = 30;

const vector<string> targetNames = {"blank", "shortCarbon", "longCarbon", "Sn112", "NatSn", "Sn124"}; 

// Extract each point from a graph and store their positions in two vectors,
// xValues and yValues
void extractGraphData(
        TGraphErrors* graph,
        vector<double>* xValues,
        vector<double>* yValues,
        vector<double>* yError)
{
    int numPoints = graph->GetN();
    
    xValues->resize(numPoints);
    yValues->resize(numPoints);
    yError->resize(numPoints);

    for(int k=0; k<numPoints; k++)
    {
        graph->GetPoint(k,xValues->at(k),yValues->at(k));
        yError->at(k) = graph->GetErrorY(k);
    }
}

bool fileExists(string fileName)
{
    ifstream infile(fileName);
    return infile.good();
}

void readGraphs(
        string runNumber,
        string driveName,
        string fileType,
        vector<string> targets,
        vector<vector<vector<double>*>*> &energies,
        vector<vector<vector<double>*>*> &crossSections,
        vector<vector<vector<double>*>*> &crossSectionsError
        )
{
    // prep vectors for filling
    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        energies.push_back(new vector<vector<double>*>);
        crossSections.push_back(new vector<vector<double>*>);
        crossSectionsError.push_back(new vector<vector<double>*>);
    }

    TFile* infile;

    // Loop over subruns to read data 
    for(int i = 0; i<=MAX_SUBRUN_NUMBER; i++)
    {
        stringstream inFileName;

        // Calculate the name of the next sub-run
        if(i < 10)
        {
            inFileName << driveName << "/analysis/run" << runNumber << "/000"
                       << i << "/" << fileType << ".root";
        }

        else if(i < 100)
        {
            inFileName << driveName << "/analysis/run" << runNumber << "/00"
                       << i << "/" << fileType << ".root";
        }

        else
        {
            cout << "Error: subrun number too large." << endl;
            exit(1);
        }

        // Attempt to open the sub-run
        if(!fileExists(inFileName.str()))
        {
            cout << "Can't open subrun " << runNumber << " " << i << endl;
            continue;
        }

        infile = new TFile(inFileName.str().c_str());
        cout << "Adding run " << runNumber << " " << i <<endl;

        // Pull out the cross section data
        for(int j=0; (size_t)j<targets.size(); j++)
        {
            TGraphErrors * graph = (TGraphErrors*)infile->Get(targets[j].c_str());
            energies[j]->push_back(new vector<double>);
            crossSections[j]->push_back(new vector<double>);
            crossSectionsError[j]->push_back(new vector<double>);
            extractGraphData(graph,energies[j]->back(),crossSections[j]->back(),crossSectionsError[j]->back());
        }

        infile->Close();
        // End of loop - move to next sub-run
    }
}

void averageGraphs(
        vector<string> targets,
        vector<vector<vector<double>*>*> &energies,
        vector<vector<vector<double>*>*> &crossSections,
        vector<vector<vector<double>*>*> &crossSectionsError
        )
{
    // create vectors for holding cross section average over all subruns:
    // crossSectionsAvg[target number]->at(data point)
    vector<vector<double>*> crossSectionsAvg;
    vector<vector<double>*> crossSectionsErrorAvg;
    vector<double> energyError;
    energyError.resize(energies[0]->at(0)->size());

    // prep vectors for filling
    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        crossSectionsAvg.push_back(new vector<double>);
        crossSectionsErrorAvg.push_back(new vector<double>);
    }

    for(int i=1; (size_t)i<energies.size(); i++)
    {
        // compute average
        crossSectionsAvg[i]->resize(energies[0]->at(0)->size());
        for(int j=0; (size_t)j<energies[1]->size(); j++)
        {
            for(int k=0; (size_t)k<energies[1]->at(0)->size(); k++)
            {
                crossSectionsAvg[i]->at(k) += crossSections[i]->at(j)->at(k);
            }
        }

        for(int k=0; (size_t)k<energies[1]->at(0)->size(); k++)
        {
            crossSectionsAvg[i]->at(k) /= energies[1]->size();
        }

        // propagate error
        crossSectionsErrorAvg[i]->resize(energies[0]->at(0)->size());
        for(int j=0; (size_t)j<energies[1]->size(); j++)
        {
            for(int k=0; (size_t)k<energies[1]->at(0)->size(); k++)
            {
                crossSectionsErrorAvg[i]->at(k) += pow(crossSectionsError[i]->at(j)->at(k),2);
            }
        }

        for(int k=0; (size_t)k<energies[1]->at(0)->size(); k++)
        {
            crossSectionsErrorAvg[i]->at(k) = pow(crossSectionsErrorAvg[i]->at(k),0.5);
            crossSectionsErrorAvg[i]->at(k) /= energies[1]->size();
        }

        // create new graphs to display the average

        TGraphErrors* graph = new TGraphErrors(energies[i]->at(0)->size(),
                                  &energies[i]->at(0)->at(0),
                                  &crossSectionsAvg[i]->at(0),
                                  &energyError[0],
                                  &crossSectionsErrorAvg[i]->at(0));
        graph->SetNameTitle(targets[i].c_str(),targets[i].c_str());
        graph->Write();
    }
}

int main(int argc, char *argv[])
{
    string runNumber = argv[1];
    string driveName = argv[2];

    // create vector for holding the energies where the cross sections
    // were calculated
    vector<vector<vector<double>*>*> energies;
    vector<vector<vector<double>*>*> energiesWaveform;

    // create vectors for holding cross section data:
    // crossSections[target number]->at(subrun number)->at(data point)
    vector<vector<vector<double>*>*> crossSections;
    vector<vector<vector<double>*>*> crossSectionsWaveform;
    vector<vector<vector<double>*>*> crossSectionsError;
    vector<vector<vector<double>*>*> crossSectionsErrorWaveform;

    readGraphs(runNumber, driveName, "cross-sections", targetNames, energies, crossSections, crossSectionsError);
    //readGraphs(runNumber, driveName, "waveform",targetNamesWaveform,energiesWaveform,crossSectionsWaveform);

    // Create output file to contain summed histos
    stringstream outfileName;
    outfileName << driveName << "/analysis/run" << runNumber << "/" << "sum.root";
    TFile *outfile = new TFile(outfileName.str().c_str(),"RECREATE");

    averageGraphs(targetNames,energies,crossSections,crossSectionsError);
    //averageGraphs(targetNamesWaveform,energiesWaveform,crossSectionsWaveform);

    outfile->Close();
}
