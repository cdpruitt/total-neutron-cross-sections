#include <iostream>
#include <string>

#include "../include/physicalConstants.h"

#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"

using namespace std;

// assuming the strongly absorbing (SA) "black body" model of the nucleus w/r/t
// incident neutrons, calculate relative difference of cross sections of two
// nuclei (mass numbers A1 and A2) at neutron energy E
double calculateCS_SA(unsigned int A, double E)
{
    double reducedWavelength = REDUCED_PLANCK_CONSTANT/pow(2*E*NEUTRON_MASS, 0.5); // in meters
    reducedWavelength *= pow(10,15)*C; // convert from units of C*s to fermi
    double nuclearRadius = NUCLEAR_RADIUS_CONSTANT*pow(A,1/(double)3); // in fermi

    return TAU*pow(nuclearRadius+reducedWavelength, 2);
}

TGraph* createSAGraph(unsigned int A1, unsigned int A2, string name)
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
}

int main()
{
    TFile* file = new TFile("ramsauer.root", "RECREATE");

    TGraph* OGraph = createSAGraph(18, 16, "oxygen");
    OGraph->SetNameTitle("relativeDiff(O18,O16)", "relativeDiff(O18,O16)");
    OGraph->Write();

    TGraph* NiGraph = createSAGraph(64, 58, "nickel");
    NiGraph->SetNameTitle("relativeDiff(Ni64,Ni58)", "relativeDiff(Ni64,Ni58)");
    NiGraph->Write();


    TGraph* SnGraph = createSAGraph(124, 112, "tin");
    SnGraph->SetNameTitle("relativeDiff(Sn124,Sn112)", "relativeDiff(Sn124,Sn112)");
    SnGraph->Write();

    file->Close();
    return 0;
}
