#ifndef OPTICAL_POTENTIAL_H
#define OPTICAL_POTENTIAL_H

#include <string>

struct OpticalPotential
{
    OpticalPotential();
    OpticalPotential(std::string configFileName);

    double VHF;
    double VHF_Edep;

    double RHF;
    double RHF_Edep;
    double RHF_Eshift;

    double aHF;

    double Avolume;
    double Avolume_Edep;
    double Avolume_Eshift;

    double Rvolume;
    double avolume;

    double Asurface;
    double Asurface_Edep;
    double Asurface_Eshift;

    double CoulombRadiusConstant;
};

#endif /* OPTICAL_POTENTIAL_H */
