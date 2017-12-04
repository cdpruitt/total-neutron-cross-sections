#include "../include/reaction.h"
#include "../include/physicalConstants.h"
#include "../include/OpticalPotential.h"

#include <iostream>

#include "TGraph.h"
#include "TFile.h"

using namespace std;

int main()
{
    vector<double> energy;
    vector<double> crossSectionValues;

    double startingPoint = log(1);
    double stepSize = (log(500)-log(1))/100;

    OpticalPotential OP("config/opticalPotential.txt");

    /*cout << "OP.VHF = " << OP.VHF << endl;
    cout << "OP.VHF_Edep = " << OP.VHF_Edep << endl;
    cout << "OP.RHF = " << OP.RHF << endl;
    cout << "OP.RHF_Edep = " << OP.RHF_Edep << endl;
    cout << "OP.RHF_Eshift = " << OP.RHF_Eshift << endl;
    cout << "OP.aHF = " << OP.aHF << endl;
    cout << "OP.Avolume = " << OP.Avolume << endl;
    cout << "OP.Avolume_Edep = " << OP.Avolume_Edep << endl;
    cout << "OP.Avolume_Eshift = " << OP.Avolume_Eshift << endl;
    cout << "OP.Rvolume = " << OP.Rvolume << endl;
    cout << "OP.avolume = " << OP.avolume << endl;
    cout << "OP.Asurface = " << OP.Asurface << endl;
    cout << "OP.Asurface_Edep = " << OP.Asurface_Edep << endl;
    cout << "OP.Asurface_Eshift = " << OP.Asurface_Eshift << endl;
    cout << "OP.CoulombRadiusConstant = " << OP.CoulombRadiusConstant << endl;
    */

    for(unsigned int i=0; i<100; i++)
    {
        double Elab = exp(i*stepSize+startingPoint);
        energy.push_back(Elab);
  
        reaction Reaction;

        Reaction.Zp = 0;
        Reaction.Z = 28;
        Reaction.A = 64;
        Reaction.DOM = 0;
        string *title = new string (" ");
        bool bool1 = 0;
        Reaction.scatter = new scat(Reaction.Zp,Reaction.Z,Reaction.A,bool1,title);
        Reaction.scatter->initIntegration();

        Reaction.Rc = OP.CoulombRadiusConstant*pow(Reaction.A,1./3.);
        Reaction.VHF = OP.VHF + OP.VHF_Edep*Elab;
        Reaction.RHF = (OP.RHF + OP.RHF_Edep*(Elab-OP.RHF_Eshift))
            *pow(Reaction.A,1./3.);
        Reaction.aHF = OP.aHF;

        Reaction.Avolume = OP.Avolume + OP.Avolume_Edep
            *(Elab-OP.Avolume_Eshift);
        Reaction.Rvolume = OP.Rvolume*pow(Reaction.A,1./3.);
        Reaction.avolume = OP.avolume;

        Reaction.Asurface = OP.Asurface + OP.Asurface_Edep*(Elab-OP.Asurface_Eshift);
        Reaction.Rsurface = Reaction.Rvolume;
        Reaction.asurface = Reaction.avolume;

        Reaction.Vso = 0.;
        Reaction.Rso = Reaction.RHF;
        Reaction.aso = Reaction.aHF;
        Reaction.AWso = 0.;

        double Ecm = Reaction.energyLab2Cm(Elab);
        Reaction.loadOM();
        Reaction.InitializeForEcm(Ecm,Elab);
        Reaction.scatter->integrateWave();

        for (int j=0;j<10;j++)
        {
            double b = ((double)j+0.5)/Reaction.scatter->Kwave;
            Reaction.scatter->getSmatrixEikonal(b);
            Reaction.scatter->TransCoef(j,(double)j+0.5);
            Reaction.scatter->getSmatrixEikonal(b);
        }

        crossSectionValues.push_back(Reaction.scatter->TotXsection()/1000);

        cout << "Lab energy = " << energy.back() << ", total cross section = " << crossSectionValues.back() << endl;
    }

    TFile* file = new TFile("~/neutronTCS/analysis/opticalModel.root","UPDATE");

    TGraph* crossSection = new TGraph(energy.size(), &energy[0], &crossSectionValues[0]);
    crossSection->SetNameTitle("Ni64","Ni64");
    crossSection->Write();

    file->Close();
}
