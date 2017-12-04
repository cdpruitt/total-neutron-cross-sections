#ifndef spinorbit_
#define spinorbit_

#include "potPara.h"
#include "imaginaryForm.h"
#include <iostream>

using namespace std;


/**
 *\brief spin orbit optical potential
 *
 * calculates the real and imaginary spin orbit potential
 */


class spinOrbit
{
 public:
  double Vzero;
  double V;
  double AW;
  double BW;
  double W;
  double DerV;
  double e;
  double R;
  double a;
  double Efermi;
  double Ecm;
  potPara Real;
  potPara Imag;
  potPara DerDispersive;
  potPara DerParticleHole;
  imaginaryForm  Form;
  void load(double,double,double,double,double,double,double);
  void load(double V0, double W0, double R0, double a0);
  void SetEnergy(double);
  double RealPotential(double,double);
  double ImaginaryPotential(double,double);
  double DerivativeDispersive(double,double);  
  double DerivativeParticleHole(double,double);
};

#endif
