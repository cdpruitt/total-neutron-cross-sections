#include "surfaceFH.h"


/**
 * loads in all the paramters defining the radial and energ dependence of the
 * surface imaginary potential
\param R0 is the radius parameter in fm
\param a0 is the diffuseness parameter in fm
\param A10 defines the magnitude of the potential
 */

void surfaceFH::load(double R0, double a0, double A0, double B0 ,
		     double C0,  double D0,  double m0, double Efermi0)

{
  R = R0;
  a = a0;
  A = A0;
  B = B0;
  C = C0;
  D = D0;
  m = m0;
  Efermi = Efermi0;

  FH.init(A,B,C,D,m,Efermi);

}
//***************************************************************
  /**
   * load parameter where the totla strength is specified not the 
   * energy dependence. THis is useful for OM fits to single energy data
   \param strength is the strength of the surface imaginary potential
   \param R0 is the radius of the potential in fm
   \param a0 is the diffuseness of the potential in fm
  */
void surfaceFH::load(double strength, double R0, double a0)
{
  R = R0;
  a = a0;
  Imaginary.init(strength,R,a);
  Dispersive.init(0.,R,a);
  DerDispersive.init(0.,R,a);
  ParticleHole.init(0.,R,a);
}
//***********************************************************
  /**
   * for each given energy this needs to be runs before calling the 
   *subsequent functions with give radial dependencies.
   */

void surfaceFH::SetEnergy(double Ecm0)
{
  Ecm = Ecm0;
  double V = FH.funct(Ecm);


  Imaginary.init(V,R,a);
  V = FH.deltaV(Ecm);
  Dispersive.init(V,R,a);
  V = FH.derDeltaV();
  DerDispersive.init(V,R,a);

  V = FH.derDeltaVParticleHole(); 
  ParticleHole.init(V,R,a);
}
//***************************************************************
/**
 * returns the surface imaginary potential at a given radius.
 * the fubction setEnergy(Ecm) must be run before using this function
\param r is radial distance in fm
 */
double surfaceFH::ImaginaryPot(double r)
{


  return Imaginary.DerWoodSaxon(r);
}
//***************************************************************
  /**
   *returns the real dispersive correction to the surface imaginary
   *potential. The function setEnergy(Ecm) must be run beforehand
   \param r is the radial distance in fm
  */
double surfaceFH::DispersiveCorrection(double r)
{
  return Dispersive.DerWoodSaxon(r);
}
//********************************************************
  /**
   *returns the derivative with respect to energy of the dispersive correction
\param r is the radial distance in fm
   */
double surfaceFH::DerivativeDispersive(double r)
{
  return DerDispersive.DerWoodSaxon(r);
}
//**********************************************************
  /**
   *returns the derivative of the hole or particle contribution to 
   *the dispersive corrections
   *for E > Efermi - particle contribution
   *for E < Efermi - hole contribution
   *It us used to calculate occupation probabilities
\param r is the radial distance in fm
   */
double surfaceFH::DerivativeParticleHole(double r)
{
  return ParticleHole.DerWoodSaxon(r);
}
//********************************************
  /**
   * returns the magnitude of the imaginary potential
     \param Ecm in the center-of-mass energy in MeV
   */
double surfaceFH::CentralImaginaryPotential(double Ecm)
{
  return FH.funct(Ecm);
}
//*******************************************
  /**
   * returns the magnitude of the dispersive correction
     \param Ecm in the center-of-mass energy in MeV
   */
double surfaceFH::CentralDeltaRealPotential(double Ecm)
{
  return FH.deltaV(Ecm);
}
//****************************************
  /**
   *returns the magnitude of the derivative of the hole 
   *or particle contribution to 
   *the dispersive corrections
   *for E > Efermi - particle contribution
   *for E < Efermi - hole contribution
   *It us used to calculate occupation probabilities
   \param Ecm is the center-of-mass energy in MeV
  */

double surfaceFH::CentralDerDeltaRealPotential(double Ecm)
{
  return FH.derDeltaV(Ecm);
}
