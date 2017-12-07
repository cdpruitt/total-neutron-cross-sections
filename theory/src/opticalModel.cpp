#include "../../analysis/include/crossSection.h"
#include "../../analysis/include/dataSet.h"
#include "../../analysis/include/dataPoint.h"
#include "../../analysis/include/CSUtilities.h"

#include "../include/physicalConstants.h"
#include "../include/OpticalPotential.h"
#include "../include/Nucleus.h"
#include "../include/calculateCS.h"

#include "TGraph.h"
#include "TFile.h"

#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char** argv)
{
    // read config files to see which nuclei should be used for calculation
    string configFileName = "config/filesToRead.txt";
    ifstream configFile(configFileName);

    if(!configFile.is_open())
    {
        cerr << "Error: failed to open " << configFileName << "; exiting..." << endl;
        exit(1);
    }

    string relativeFileName = "config/relativePlots.txt";
    ifstream relativeFile(relativeFileName);

    if(!relativeFile.is_open())
    {
        cerr << "Error: failed to open " << relativeFileName << "; exiting..." << endl;
        exit(1);
    }

    string buffer;
    vector<string> nuclides;
    vector<pair<string,string>> relativeDiffNames;

    while(getline(configFile, buffer))
    {
        nuclides.push_back(buffer);
    }

    while(getline(relativeFile, buffer))
    {
        vector<string> tokens;
        istringstream iss(buffer);
        copy(istream_iterator<string>(iss),
                istream_iterator<string>(),
                back_inserter(tokens));

        relativeDiffNames.push_back(make_pair(tokens[0], tokens[1]));
    }

    // begin cross section calculation
    vector<CrossSection> crossSections;

    for(auto& nucleus : nuclides)
    {
        Nucleus N("config/" + nucleus + ".txt");

        cout << "Producing cross section for " << N.Name << endl;

        // define energy range
        vector<double> energyRange;

        double startingPoint = log(1);
        double stepSize = (log(500)-log(1))/100;

        for(unsigned int i=0; i<100; i++)
        {
            energyRange.push_back(exp(i*stepSize+startingPoint));
        }

        // identify which optical model should be used for cross section
        // calculation, and perform calculation
        string argument = argv[1];

        if(argument == "SA")
        {
            // strongly absorbing disk model

            crossSections.push_back(calculateCS_SA(N, energyRange));
        }

        else if(argument == "ROP")
        {
            // optical potential with only a real, Woods-Saxon part

            OpticalPotential OP("config/" + nucleus + "_ROP.txt");
            crossSections.push_back(calculateCS_ROP(N, OP, energyRange));
        }

        else if(argument == "COP")
        {
            // optical potential with real and imaginary Woods-Saxon parts

            OpticalPotential OP("config/" + nucleus + "_COP.txt");
            crossSections.push_back(calculateCS_COP(N, OP, energyRange));
        }

        else if(argument == "COPS")
        {
            // optical potential with real and imaginary Woods-Saxon and surface
            // (derivative of Woods-Saxon) parts

            OpticalPotential OP("config/" + nucleus + "_COPS.txt");
            crossSections.push_back(calculateCS_COPS(N, OP, energyRange));
        }

        else
        {
            cerr << "Error: argument describing the optical potential to be used must be supplied." << endl;
        return 1;
        }
    }

    // cross sections have been calculated; make plots
    TFile* file = new TFile("opticalModel.root","RECREATE");

    for(auto& cs : crossSections)
    {
        TGraph* graph = new TGraph(
                cs.getDataSet().getNumberOfPoints(),
                &cs.getDataSet().getXValues()[0],
                &cs.getDataSet().getYValues()[0]);

        graph->SetNameTitle(cs.name.c_str(), cs.name.c_str());
        graph->Write();
    }

    for(auto& pair : relativeDiffNames)
    {
        CrossSection firstCS;
        CrossSection secondCS;

        for(auto& cs : crossSections)
        {
            if(cs.name == pair.first)
            {
                firstCS = cs;
            }

            if(cs.name == pair.second)
            {
                secondCS = cs;
            }
        }

        if(firstCS.name=="" || secondCS.name=="")
        {
            cerr << "Error: failed to find both cross sections \""
                << pair.first << "\", \"" << pair.second
                << " needed to calculate relative difference." << endl;

            break;
        }

        CrossSection relativeDiff = calculateRelativeDiff(firstCS, secondCS);
        relativeDiff.name = "relativeDiff(" + pair.first + "," + pair.second + ")";
        TGraph* graph = new TGraph(
                relativeDiff.getDataSet().getNumberOfPoints(),
                &relativeDiff.getDataSet().getXValues()[0],
                &relativeDiff.getDataSet().getYValues()[0]);

        graph->SetNameTitle(relativeDiff.name.c_str(), relativeDiff.name.c_str());
        graph->Write();
    }

    file->Close();

    return 0;
}
