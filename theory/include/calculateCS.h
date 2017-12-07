#ifndef CALCULATE_CS_H
#define CALCULATE_CS_H

#include "../../analysis/include/crossSection.h"
#include "Nucleus.h"
#include "OpticalPotential.h"

#include <vector>

CrossSection calculateCS_SA(
        Nucleus N,
        std::vector<double> energyRange
        );

CrossSection calculateCS_ROP(
        Nucleus N,
        OpticalPotential OP,
        std::vector<double> energyRange
        );

CrossSection calculateCS_COP(
        Nucleus N,
        OpticalPotential OP,
        std::vector<double> energyRange
        );

CrossSection calculateCS_COPS(
        Nucleus N,
        OpticalPotential OP,
        std::vector<double> energyRange
        );


#endif /* CALCULATE_CS_H */
