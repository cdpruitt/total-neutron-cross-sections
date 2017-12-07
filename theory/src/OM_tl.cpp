#include "../include/reaction.h"
#include <iostream>

using namespace std;

CrossSection calculateCS(
        Nucleus N,
        OpticalPotential OP,
        vector<double> energyRange
        )
{
    reaction Reaction;

    Reaction.Zp = 0;
    Reaction.Z = N.Z;
    Reaction.A = N.A;
    Reaction.DOM = 0;
    string *title = new string (" ");

    Reaction.scatter = new scat(Reaction.Zp,Reaction.Z,Reaction.A,false,title);
    Reaction.scatter->initIntegration();

    DataSet dataSet;

    for(auto& energy : energyRange) // energy in the lab frame expected
    {
        double energyCM = Reaction.energyLab2Cm(energy);

        Reaction.Rc = OP.CoulombRadiusConstant*pow(Reaction.A,1./3.);
        Reaction.VHF = OP.VHF + OP.VHF_Edep*energy;
        Reaction.RHF = (OP.RHF + OP.RHF_Edep*(energy-OP.RHF_Eshift))
            *pow(Reaction.A,1./3.);
        Reaction.aHF = OP.aHF;

        Reaction.Avolume = OP.Avolume + OP.Avolume_Edep
            *(energy-OP.Avolume_Eshift);
        Reaction.Rvolume = OP.Rvolume*pow(Reaction.A,1./3.);
        Reaction.avolume = OP.avolume;

        Reaction.Asurface = OP.Asurface + OP.Asurface_Edep*(energy-OP.Asurface_Eshift);
        Reaction.Rsurface = Reaction.Rvolume;
        Reaction.asurface = Reaction.avolume;

        Reaction.Vso = 0.;
        Reaction.Rso = Reaction.RHF;
        Reaction.aso = Reaction.aHF;
        Reaction.AWso = 0.;

        // loop over energies
        Reaction.loadOM();
        Reaction.InitializeForEcm(energyCM, energy);
        Reaction.scatter->integrateWave();

        for (int j=0;j<10;j++)
        {
            double b = ((double)j+0.5)/Reaction.scatter->Kwave;
            Reaction.scatter->getSmatrixEikonal(b);
            Reaction.scatter->TransCoef(j,(double)j+0.5);
            Reaction.scatter->getSmatrixEikonal(b);
        }

        double value = Reaction.scatter->TotXsection()/1000;

        DataPoint dataPoint;

        dataPoint.setXValue(energy);
        dataPoint.setYValue(value);

        dataSet.addPoint(dataPoint);
    }

    CrossSection cs(N.Name);
    cs.addDataSet(dataSet);

    return cs;
}
