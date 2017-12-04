#ifndef surfaceTF_
#define surfaceTF_

#include "potPara.h"
#include "twoFermi.h"
#include <iostream>

using namespace std;

/**
 *\brief imaginary surface potential and dispersive correction
 *
 *this class deals with all aspects of the surface imaginary potential, its
 *dispersive correction and contributions to effective mass and occupation
 *probabilities. the magnitude of the imaginary potential is parametrized as
 *\f$ W(E) = A \frac{\|X_{1}\|^{m}}{\|X_{1}\|^{m}+B^m} \frac{1}{1+\exp\left( \frac{X_{1}-C}{D}\right)} \frac{1}{1+\exp-\left(\frac{X_{1}+C}{D}\right)} \f$ 
 * \f$ X_{1} = E-E_{Fermi} \f$
 */


class surfaceTF
{
 public:
  void load(double,double,double,double,double,double,double,double);
  void load(double strength, double R0, double a0);
  void SetEnergy(double);
  double ImaginaryPot(double);
  double DispersiveCorrection(double);
  double DerivativeDispersive(double);
  double DerivativeParticleHole(double);

  double CentralImaginaryPotential(double);
  double CentralDeltaRealPotential(double);
  double CentralDerDeltaRealPotential(double);
  double CentralParticleHole(double);
  double getMaxW();

  potPara Imaginary;  //imaginary potential 
  potPara Dispersive; //dispersive coorection to real potential
  potPara DerDispersive; //derivative of dispersive need for effective mass
  potPara ParticleHole; // derivative of particle of hole contribution to 
                   // dispersive correction needed for occupation probabilities

  twoFermi TF; // energy dependence of imaginary potential

  double A;
  double B;
  double C;
  double D;
  double Wstart;
  double Efermi;
  double Ecm;

  double R;
  double a;

};

#endif 
