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

// assuming the strongly absorbing "black body" model of the nucleus w/r/t
// incident neutrons, calculate relative difference
double calculateRelDiff_SA(unsigned int A1, unsigned int A2, double E)
{
    return (calculateCS_SA(A1,E)-calculateCS_SA(A2,E))/(calculateCS_SA(A1,E)+calculateCS_SA(A2,E));
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
    OGraph->SetLineColor(kBlue);
    OGraph->SetLineWidth(2);
    OGraph->SetLineStyle(1);

    TGraph* NiGraph = createSAGraph(64, 58, "nickel");
    NiGraph->SetLineColor(kGreen);
    NiGraph->SetLineWidth(2);
    NiGraph->SetLineStyle(1);

    TGraph* SnGraph = createSAGraph(124, 112, "tin");
    SnGraph->SetLineColor(kRed);
    SnGraph->SetLineWidth(2);
    SnGraph->SetLineStyle(1);

    TMultiGraph* allGraphs = new TMultiGraph();
    allGraphs->SetNameTitle("stronglyAbsorbing","stronglyAbsorbing");

    allGraphs->Add(OGraph,"AL");
    allGraphs->Add(NiGraph,"AL");
    allGraphs->Add(SnGraph,"AL");
    allGraphs->Write();

    allGraphs->Draw("AL");

    file->Close();
    return 0;
}
