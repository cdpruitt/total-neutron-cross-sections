#ifndef PHYSICAL_CONSTANTS_H
#define PHYSICAL_CONSTANTS_H

#include <math.h>

/******************************************************************************/
/* Physical constants used in the experiment */
/******************************************************************************/

const double C = 299792458; // speed of light in m/s
const double NEUTRON_MASS = 939.56541; // in MeV/c^2
const double AVOGADROS_NUMBER = 6.02214*pow(10.,23.); // in atoms/mol

const double REDUCED_PLANCK_CONSTANT = 6.5821195*pow(10,-22); // in MeV*s
const double NUCLEAR_RADIUS_CONSTANT = 1.07; // in fermi, for R = R0 * A^(1/3)

const double TAU = 6.2831853; // = 2*pi

#endif
