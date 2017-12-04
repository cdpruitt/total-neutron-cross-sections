#ifndef hartreefock_
#define hartreefock_

#include "potPara.h"
#include <cmath>
#include <iostream>

using namespace std;

/**
 *\brief give the energy-dependent Hartree-Fock potential
 *
 *Deals with all aspects of the Hartree-Fock potential - This is the
 * main real potential without the dispersive corrections. The energy
 * dependece is to take into account nonlocality
 */
class hartreeFock
{
public:

  void load(double,double,double,double,double,double,double,
            double,double,double);
  void load(double V0, double dV0, double R0, double a0);
  void SetEnergy(double);
  double RealPotential(double);
  double DerivativePotential(double);

  potPara RealVolume;
  potPara DerivativeRealVolume;
  potPara surface;
  potPara DerivativeSurface;
  double RealSurface(double Ecm);

  double R;
  double a;
  double VHF;
  double alpha;
  double beta;
  double gamma;
  double Efermi;
  double alphaS;
  double betaS;
  double gammaS;

  double V; //!< depth of hartree fock potential
  double dV; //!< derivarive of V with respect to com energy

  static double const dd;
  static double const ss;

};

#endif
