#include <iostream>
#include <string>

#include "../include/physicalConstants.h"
#include "../include/ramsauer.h"

#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"

using namespace std;

// assuming the strongly absorbing (SA) "black body" model of the nucleus w/r/t
// incident neutrons, calculate cross sections of a nucleus (mass numbers A) at neutron energy E
double calculateCS_SA(unsigned int A, double E)
{
    double reducedWavelength = (REDUCED_PLANCK_CONSTANT*((A+1)/A))/pow(2*E*NEUTRON_MASS, 0.5); // in meters
    reducedWavelength *= pow(10,15)*C; // convert from units of C*s to fermi
    double nuclearRadius = NUCLEAR_RADIUS_CONSTANT*pow(A,1/(double)3); // in fermi

    double crossSection = TAU*pow(nuclearRadius+reducedWavelength, 2); // in fermi^2
    return crossSection*FM2_TO_BARNS; // in barns
}

TGraph* createSAGraph(double A, string name)
{
    vector<double> energy;
    vector<double> crossSection;

    double startingPoint = log(1);
    double stepSize = (log(500)-log(1))/100;

    for(unsigned int i=0; i<100; i++)
    {
        double E = exp(i*stepSize+startingPoint);
        energy.push_back(E);
        crossSection.push_back(calculateCS_SA(A, E));
    }

    TGraph* SAGraph = new TGraph(energy.size(), &energy[0], &crossSection[0]);
    SAGraph->SetNameTitle(name.c_str(),name.c_str());
    return SAGraph;
}

/*TGraph* createSARelDiffGraph(double A1, double A2, string name)
{
    vector<double> energy;
    vector<double> relDiff;

    double startingPoint = log(1);
    double stepSize = (log(500)-log(1))/100;

    for(unsigned int i=0; i<100; i++)
    {
        double E = exp(i*stepSize+startingPoint);
        energy.push_back(E);
        relDiff.push_back(calculateRelDiff_SA(A1, A2, E));
    }

    TGraph* relDiff_SA = new TGraph(energy.size(), &energy[0], &relDiff[0]);
    relDiff_SA->SetNameTitle(name.c_str(),name.c_str());
    return relDiff_SA;
}*/

int main()
{
    TFile* file = new TFile("ramsauer.root", "RECREATE");

    TGraph* CGraph = createSAGraph(12, "carbon");
    CGraph->SetNameTitle("SA_A=12", "SA_A=12");
    CGraph->Write();

    TGraph* NiGraph = createSAGraph(58.7, "nickel");
    NiGraph->SetNameTitle("SA_A=58.7", "SA_S=58.7");
    NiGraph->Write();

    TGraph* SnGraph = createSAGraph(118.7, "tin");
    SnGraph->SetNameTitle("SA_A=118.7", "SA_A=118.7");
    SnGraph->Write();

    TGraph* PbGraph = createSAGraph(207.2, "lead");
    PbGraph->SetNameTitle("SA_A=207.2", "SA_A=207.2");
    PbGraph->Write();

    file->Close();
    return 0;
}
