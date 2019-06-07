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
    double relativisticP = pow(pow(E/C,2) + 2*E*NEUTRON_MASS, 0.5); // in MeV/C
    double reducedWavelength = (REDUCED_PLANCK_CONSTANT*((A+1)/A))/relativisticP; // in meters
    reducedWavelength *= pow(10,15)*C; // convert from units of C*s to fermi
    double nuclearRadius = NUCLEAR_RADIUS_CONSTANT*pow(A,1/(double)3); // in fermi

    double crossSection = TAU*pow(nuclearRadius+reducedWavelength, 2); // in fermi^2
    return crossSection*FM2_TO_BARNS; // in barns
}

double calculateCS_Ramsauer(unsigned int A, double E)
{
    double relativisticP = pow(pow(E/C,2) + 2*E*NEUTRON_MASS, 0.5); // in MeV/C
    double reducedWavelength = (REDUCED_PLANCK_CONSTANT*((A+1)/A))/relativisticP; // in meters
    reducedWavelength *= pow(10,15)*C; // convert from units of C*s to fermi
    //cout << "For E = " << E << ", reduced wavelength = " << reducedWavelength << endl;

    double nuclearRadius = NUCLEAR_RADIUS_CONSTANT*pow(A,1/3.0); // in fermi
    double delta = 4/3.0 * nuclearRadius * (pow((E+WOODS_SAXON_WELL_DEPTH)/E,0.5)-1)/reducedWavelength;
    double deltaStar = 2 * SKIN_THICKNESS * (pow((E+WOODS_SAXON_WELL_DEPTH)/E,0.5)-pow((E+WOODS_SAXON_WELL_DEPTH/2)/E,0.5))/reducedWavelength;
    //double deltaStar = 1.25;
    //double delta = 2.18*pow(A,1/3.0) * (pow((E+WOODS_SAXON_WELL_DEPTH)/E,0.5)-1)/reducedWavelength;

    double crossSection = TAU*pow(nuclearRadius+reducedWavelength, 2)*(1-0.104*cos(delta-deltaStar)); // in fermi^2
    return crossSection*FM2_TO_BARNS; // in barns
}

double calculateCS_OM(unsigned int A, double E)
{
    double reducedWavelength = (REDUCED_PLANCK_CONSTANT*((A+1)/A))/pow(2*E*NEUTRON_MASS, 0.5); // in meters
    reducedWavelength *= pow(10,15)*C; // convert from units of C*s to fermi

    double nuclearRadius = NUCLEAR_RADIUS_CONSTANT*pow(A,1/(double)3); // in fermi

    const double U = 200; // in MeV
    double phaseShift = (4*nuclearRadius/3.)*(pow((E+U)/E,0.5)-1)/reducedWavelength; // unitless

    double phaseOffset = -0.4;

    const double rho = 0.15; // imaginary component

    double crossSection = TAU*pow(nuclearRadius+reducedWavelength,2)
        *(1-rho*cos(phaseShift+phaseOffset));
    return crossSection*FM2_TO_BARNS; // in barns
}

TGraph* createSAGraph(double A, string name)
{
    vector<double> energy;
    vector<double> crossSection;

    double startingPoint = log(1);
    double stepSize = (log(600)-log(1))/100;

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

TGraph* createRamsauerGraph(double A, string name)
{
    vector<double> energy;
    vector<double> crossSection;

    double startingPoint = log(1);
    double stepSize = (log(600)-log(1))/100;

    for(unsigned int i=0; i<100; i++)
    {
        double E = exp(i*stepSize+startingPoint);
        energy.push_back(E);
        crossSection.push_back(calculateCS_Ramsauer(A, E));
    }

    TGraph* RamsauerGraph = new TGraph(energy.size(), &energy[0], &crossSection[0]);
    RamsauerGraph->SetNameTitle(name.c_str(),name.c_str());
    return RamsauerGraph;
}

TGraph* createOMGraph(double A, string name)
{
    vector<double> energy;
    vector<double> crossSection;

    double startingPoint = log(1);
    double stepSize = (log(600)-log(1))/100;

    for(unsigned int i=0; i<100; i++)
    {
        double E = exp(i*stepSize+startingPoint);
        energy.push_back(E);
        crossSection.push_back(calculateCS_OM(A, E));
    }

    TGraph* OMGraph = new TGraph(energy.size(), &energy[0], &crossSection[0]);
    OMGraph->SetNameTitle(name.c_str(),name.c_str());
    return OMGraph;
}

TGraph* createRelDiffGraph(TGraph* graph1, TGraph* graph2, string name)
{
    vector<double> energy1;
    vector<double> energy2;

    vector<double> crossSection1;
    vector<double> crossSection2;

    for(int i=1; i<graph1->GetN(); i++)
    {
        double xValue = 0;
        double yValue = 0;

        graph1->GetPoint(i, xValue, yValue);

        energy1.push_back(xValue);
        crossSection1.push_back(yValue);
    }

    for(int i=1; i<graph2->GetN(); i++)
    {
        double xValue = 0;
        double yValue = 0;

        graph2->GetPoint(i, xValue, yValue);

        energy2.push_back(xValue);
        crossSection2.push_back(yValue);
    }

    vector<double> relDiff;

    for(int i=0; i<crossSection1.size(); i++)
    {
        double RD = (crossSection1[i]-crossSection2[i])/(crossSection1[i]+crossSection2[i]);
        relDiff.push_back(100*RD); // in percent
    }

    TGraph* relDiffGraph = new TGraph(energy1.size(), &energy1[0], &relDiff[0]);

    return relDiffGraph;
}

int main()
{
    TFile* file = new TFile("ramsauer.root", "RECREATE");

    TGraph* CGraph = createSAGraph(12, "carbon");
    CGraph->SetNameTitle("SA_A=12", "SA_A=12");
    CGraph->Write();

    TGraph* O16Graph = createSAGraph(16, "oxygen16");
    O16Graph->SetNameTitle("SA_A=16", "SA_S=16");
    O16Graph->Write();

    TGraph* O18Graph = createSAGraph(18, "oxygen18");
    O18Graph->SetNameTitle("SA_A=18", "SA_S=18");
    O18Graph->Write();

    TGraph* NiGraph = createSAGraph(58.7, "nickel");
    NiGraph->SetNameTitle("SA_A=58.7", "SA_S=58.7");
    NiGraph->Write();

    TGraph* Ni58Graph = createSAGraph(58, "nickel58");
    Ni58Graph->SetNameTitle("SA_A=58", "SA_S=58");
    Ni58Graph->Write();

    TGraph* Ni64Graph = createSAGraph(64, "nickel64");
    Ni64Graph->SetNameTitle("SA_A=64", "SA_S=64");
    Ni64Graph->Write();

    TGraph* SnGraph = createSAGraph(118.7, "tin");
    SnGraph->SetNameTitle("SA_A=118.7", "SA_A=118.7");
    SnGraph->Write();

    TGraph* Sn112Graph = createSAGraph(112, "tin112");
    Sn112Graph->SetNameTitle("SA_A=112", "SA_S=112");
    Sn112Graph->Write();

    TGraph* Sn124Graph = createSAGraph(124, "tin124");
    Sn124Graph->SetNameTitle("SA_A=124", "SA_S=124");
    Sn124Graph->Write();

    TGraph* SnGraphOM = createOMGraph(118.7, "tin, OM");
    SnGraphOM->SetNameTitle("OM_A=118.7", "OM_A=118.7");
    SnGraphOM->Write();

    TGraph* PbGraph = createSAGraph(207.2, "lead");
    PbGraph->SetNameTitle("SA_A=207.2", "SA_A=207.2");
    PbGraph->Write();

    TGraph* CGraphRamsauer = createRamsauerGraph(12, "carbon");
    CGraphRamsauer->SetNameTitle("Ramsauer_A=12", "Ramsauer_A=12");
    CGraphRamsauer->Write();

    //file->Close();
    //exit(0);

    TGraph* O16GraphRamsauer = createRamsauerGraph(16, "oxygen16");
    O16GraphRamsauer->SetNameTitle("Ramsauer_A=16", "Ramsauer_S=16");
    O16GraphRamsauer->Write();

    TGraph* O18GraphRamsauer = createRamsauerGraph(18, "oxygen18");
    O18GraphRamsauer->SetNameTitle("Ramsauer_A=18", "Ramsauer_S=18");
    O18GraphRamsauer->Write();

    TGraph* NiGraphRamsauer = createRamsauerGraph(58.7, "nickel");
    NiGraphRamsauer->SetNameTitle("Ramsauer_A=58.7", "Ramsauer_S=58.7");
    NiGraphRamsauer->Write();

    TGraph* Ni58GraphRamsauer = createRamsauerGraph(58, "nickel58");
    Ni58GraphRamsauer->SetNameTitle("Ramsauer_A=58", "Ramsauer_S=58");
    Ni58GraphRamsauer->Write();

    TGraph* Ni64GraphRamsauer = createRamsauerGraph(64, "nickel64");
    Ni64GraphRamsauer->SetNameTitle("Ramsauer_A=64", "Ramsauer_S=64");
    Ni64GraphRamsauer->Write();

    TGraph* SnGraphRamsauer = createRamsauerGraph(118.7, "tin");
    SnGraphRamsauer->SetNameTitle("Ramsauer_A=118.7", "Ramsauer_A=118.7");
    SnGraphRamsauer->Write();

    TGraph* Sn112GraphRamsauer = createRamsauerGraph(112, "tin112");
    Sn112GraphRamsauer->SetNameTitle("Ramsauer_A=112", "Ramsauer_S=112");
    Sn112GraphRamsauer->Write();

    TGraph* Sn124GraphRamsauer = createRamsauerGraph(124, "tin124");
    Sn124GraphRamsauer->SetNameTitle("Ramsauer_A=124", "Ramsauer_S=124");
    Sn124GraphRamsauer->Write();

    TGraph* PbGraphRamsauer = createRamsauerGraph(207.2, "lead");
    PbGraphRamsauer->SetNameTitle("Ramsauer_A=207.2", "Ramsauer_A=207.2");
    PbGraphRamsauer->Write();

    TGraph* ORelDiffGraph = createRelDiffGraph(createSAGraph(18, "O18"), createSAGraph(16, "O16"), "ORelDiff");
    ORelDiffGraph->SetNameTitle("RelDiff18_16", "RelDiff18_16");
    ORelDiffGraph->Write();

    TGraph* NiRelDiffGraph = createRelDiffGraph(createSAGraph(64, "Ni64"), createSAGraph(58, "Ni58"), "NiRelDiff");
    NiRelDiffGraph->SetNameTitle("RelDiff64_58", "RelDiff64_58");
    NiRelDiffGraph->Write();

    TGraph* SnRelDiffGraphThird = createRelDiffGraph(createSAGraph(124, "Sn124"), createSAGraph(112, "Sn112"), "SnRelDiff");
    SnRelDiffGraphThird->SetNameTitle("RelDiff124_112Third", "RelDiff124_112Third");
    SnRelDiffGraphThird->Write();

    TGraph* SnRelDiffGraphSixth = createRelDiffGraph(createSAGraph(118, "Sn124"), createSAGraph(112, "Sn112"), "SnRelDiff");
    SnRelDiffGraphSixth->SetNameTitle("RelDiff124_112Sixth", "RelDiff124_112Sixth");
    SnRelDiffGraphSixth->Write();

    TGraph* ORelDiffRamsauerGraph = createRelDiffGraph(createRamsauerGraph(18, "O18"), createRamsauerGraph(16, "O16"), "ORelDiffRamsauer");
    ORelDiffRamsauerGraph->SetNameTitle("RelDiffRamsauer18_16", "RelDiffRamsauer18_16");
    ORelDiffRamsauerGraph->Write();

    TGraph* NiRelDiffRamsauerGraph = createRelDiffGraph(createRamsauerGraph(64, "Ni64"), createRamsauerGraph(58, "Ni58"), "NiRelDiffRamsauer");
    NiRelDiffRamsauerGraph->SetNameTitle("RelDiffRamsauer64_58", "RelDiffRamsauer64_58");
    NiRelDiffRamsauerGraph->Write();

    TGraph* SnRelDiffRamsauerGraph = createRelDiffGraph(createRamsauerGraph(124, "Sn124"), createRamsauerGraph(112, "Sn112"), "SnRelDiffRamsauer");
    SnRelDiffRamsauerGraph->SetNameTitle("RelDiffRamsauer124_112", "RelDiffRamsauer124_112");
    SnRelDiffRamsauerGraph->Write();

    TGraph* SnRelDiffRamsauerSixthGraph = createRelDiffGraph(createRamsauerGraph(118, "Sn124Sixth"), createRamsauerGraph(112, "Sn112"), "SnRelDiffRamsauer");
    SnRelDiffRamsauerSixthGraph->SetNameTitle("RelDiffRamsauerSixth124_112", "RelDiffRamsauerSixth124_112");
    SnRelDiffRamsauerSixthGraph->Write();

    file->Close();
    return 0;
}
