#include "../include/twoFermi.h"



/**
 * Constructor
 */
twoFermi::twoFermi(double A0, double B0, double C0, double D0,
    double Wstart0,double Ef0) : disperse()
{
  init(A0,B0,C0,D0,Wstart0,Ef0);
}
//*****************************************************************
/**
 * initialization of the the class if one wants to change the parameters
 */
void twoFermi::init(double A0, double B0, double C0, double D0,
     double Wstart0, double Ef0)
{
  A = A0;
  B = B0;
  C = C0;
  D = D0;
  Wstart = Wstart0;
  Ef = Ef0;
  halfIntegral =0;
  halfIntegral = deltaVXAboveBelow(0.);
}
//*************************************************************
  /**
   * returns the magnitude of the surface imaginary potential
   \param E is the center-of-mass energy on the nucleon in MeV
  */
double twoFermi::funct(double E)
{
  double x = E - Ef;
  return functX(x);
}
//*************************************************************
  /**
   * returns the magnitude of the surface imaginary potential
   \param x = Ecm-Efermi  in MeV
  */
double twoFermi::functX(double x)
{
  double fact = abs(x) - Wstart;
  if (fact <= 0.) return 0.;
  double fact1 = exp(fact/B);

  return A/(1.+exp((abs(x)-C)/D))*(fact1-1.)/(fact1+1.);

}
//***************************************************************
  /**
   * returns the energy derivative of the magnitude of 
   * the surface imaginary potential
   \param E in the center-of-mass energy of the nucleon in MeV
  */
double twoFermi::derFunct(double E)
{
  double x = E - Ef;
  return derFunctX(x);
}
//***************************************************************
  /**
   * returns the energy derivative of the magnitude of 
   * the surface imaginary potential
   \param x = Ecm-Efermi in MeV
  */
double twoFermi::derFunctX(double x)
{
  double fact = abs(x) - Wstart;
  if (fact <= 0.) return 0.;
  double fact1 = exp(fact/B);
  double dFact1 = fact1/B;
  double fact2 = exp((abs(x)-C)/D);
  double dFact2 = fact2/D; 

  double term1 = 1./(1.+fact2);
  double dTerm1 = -pow(term1,2)*dFact2;
  double term2 = fact1 - 1.;
  double dTerm2 = dFact1;
  double term3 = 1./(1.+fact1);
  double dTerm3 = -pow(term3,2)*dFact1;

  double out = A*(dTerm1*term2*term3 + term1*dTerm2*term3
	    + term1*term2*dTerm3);
  if (x < 0.) out *= -1.;
  return out;
}
//****************************************************************
  /**
   * returns the dispersive correction
   \param E is the center-of-mass energy of the nucleon in MeV
   */
double twoFermi::deltaV(double E)
{
  return deltaVX(E-Ef);
}
//*************************************************************
  /**
   * returns the derivative of the real dispersive correction associated with
   * the surface imaginary potential
   \param E is the center-of-mass energy of the nucleon in MeV
   */
double twoFermi::derDeltaV(double E)
{
  return derDeltaVX(E-Ef);
}
//*************************************************************
  /**
   * returns the derivative of the real dispersive correction associated with
   * the surface imaginary potential. The function deltaV must be run first.
   */
double twoFermi::derDeltaV()
{
  return derDeltaVX();
}

//*************************************************************
  /**
   * returns the disperive factor used for Hole occupation probabilities
   \param E is the center-of-mass energy in MeV
  */
double twoFermi::derDeltaVHole(double E)
{
  return derDeltaVHoleX(E-Ef);
}
//*************************************************************
  /**
   * returns the disperive factor used for particle occupation probabilities
   \param E is the center-of-mass energy in MeV
  */
double twoFermi::derDeltaVParticle(double E)
{
  return derDeltaVParticleX(E-Ef);
}
//*************************************************************
  /**
   * returns the disperive factor used for particle (hole)occupation 
   * probabilities when the energy is above (below) the fermi energy.
   * The function deltaV must be run first.
  */
double twoFermi::derDeltaVParticleHole()
{
  return derDeltaVParticleHoleX();
}
//************************************************
double twoFermi::deltaVAboveBelow(double E)
{
  return deltaVXAboveBelow(E-Ef);
}
