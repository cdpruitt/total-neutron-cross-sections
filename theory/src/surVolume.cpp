#include "../include/surVolume.h"



/**
 * constructor
 */
surVolume::surVolume(double A0, double B0, double sigma0, double m0, 
		 double Ef0, double Ea0, double alpha0,double Ep0) : disperse()
{
  init(A0,B0,sigma0,m0,Ef0,Ea0,alpha0,Ep0);
}
//*****************************************************************
/**
 * initialization of the the class if one wants to change the parameters
 */
void surVolume::init(double A0, double B0, double sigma0, double m0, 
		     double Ef0, double Ea0, double alpha0, double Ep0)
{
  A = A0;
  B = B0;
  sigma = sigma0;
  m = m0;
  Ef = Ef0;
  Ea = Ea0;
  masy = 2.;
  alpha = alpha0;
  Ep = Ep0;

}
//*************************************************************
  /**
   * returns the magnitude of the surface imaginary potential
   \param E is the center-of-mass energy on the nucleon in MeV
  */
double surVolume::funct(double E)
{
  double x = E - Ef;
  return functX(x);
}
//*************************************************************
  /**
   * returns the magnitude of the surface imaginary potential
   \param x = Ecm-Efermi  in MeV
  */
double surVolume::functX(double x)
{
  double out;
  if (abs(x) <= Ep) return 0.;
  out= A*pow(abs(x)-Ep,m)/(pow(abs(x)-Ep,m)+pow(B,m));
  if (x < -Ea) out -= A*pow(x+Ea,masy)/(pow(x+Ea,masy)+pow(Ea,masy));
    else if (x > Ea) out += alpha*(sqrt(x+Ef)+pow(Ef+Ea,3./2.)/2./(x+Ef)
				   -1.5*sqrt(Ef+Ea));

  //out *= exp(-pow(x/sigma,2)/2.);
  out *= exp(-abs(x)/sigma);
  return out; 
}
//***************************************************************
  /**
   * returns the energy derivative of the magnitude of 
   * the surface imaginary potential
   \param E in the center-of-mass energy of the nucleon in MeV
  */
double surVolume::derFunct(double E)
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
double surVolume::derFunctX(double x)
{
  double denom = pow(abs(x)-Ep,m)+pow(B,m);
  double one = pow(abs(x)-Ep,m)/denom;
  double derOne;
  if (abs(x) <= Ep) return 0.;
  else 
    {
     derOne = m/(abs(x)-Ep)*(pow(abs(x)-Ep,m)/denom
		  - pow(abs(x)-Ep,2*m)/pow(denom,2));
     if (x < 0.) derOne *= -1.;
    } 
 

  if (x < -Ea)
    {
      double denom2 = pow(x+Ea,masy) + pow(Ea,masy);
      double one2 = pow(x+Ea,masy)/denom2;
      double derOne2 =  masy/(x+Ea)*(pow(abs(x+Ea),masy)/denom2
		  - pow(abs(x+Ea),2.*masy)/pow(denom2,2));
      one -= one2;
      derOne -= derOne2;
    }
  else if (x > Ea)
    {
      one += alpha*(sqrt(x+Ef)+pow(Ef+Ea,3./2.)/2./(x+Ef)
		    -1.5*sqrt(Ef+Ea));
       derOne += alpha*(0.5/sqrt(x+Ef) - 
			pow(Ef+Ea,3./2.)/2./pow(x+Ef,2));
    }


  double two = exp(-abs(x)/sigma);
  double derTwo = -two/sigma;
  if (x < 0.) derTwo *= -1.;

  //double two = exp(-pow(x/sigma,2)/2.);
  //double derTwo = -x/pow(sigma,2)*two;


  return A*(derOne*two+one*derTwo);
}
//****************************************************************
  /**
   * returns the magnitude of the real dispersive correction associated with
   * the surface imaginary potential
   \param E is the center-of-mass energy of the nucleon in MeV
   */
double surVolume::deltaV(double E)
{
  return deltaVX(E-Ef);
}
//*************************************************************
  /**
   * returns the derivative of the real dispersive correction associated with
   * the surface imaginary potential
   \param E is the center-of-mass energy of the nucleon in MeV
   */
double surVolume::derDeltaV(double E)
{
  return derDeltaVX(E-Ef);
}
//*************************************************************
  /**
   * returns the disperive factor used for Hole occupation probabilities
   \param E is the center-of-mass energy in MeV
  */
double surVolume::derDeltaVHole(double E)
{
  return derDeltaVHoleX(E-Ef);
}
//*************************************************************
  /**
   * returns the disperive factor used for Particle occupation probabilities
   \param E is the center-of-mass energy in MeV
  */
double surVolume::derDeltaVParticle(double E)
{
  return derDeltaVParticleX(E-Ef);
}
