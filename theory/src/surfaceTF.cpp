#include "../include/surfaceTF.h"


/**
 * loads in all the paramters defining the radial and energ dependence of the
 * surface imaginary potential
\param R0 is the radius parameter in fm
\param a0 is the diffuseness parameter in fm
\param A10 defines the magnitude of the potential
 */

void surfaceTF::load(double R0, double a0, double A0, double B0 ,
		     double C0,  double D0, 
                     double Wstart0, double Efermi0)

{
  R = R0;
  a = a0;
  A = A0;
  B = B0;
  C = C0;
  D = D0;
  Wstart = Wstart0;
  Efermi = Efermi0;

  TF.init(A,B,C,D,Wstart,Efermi);

}
//***************************************************************
  /**
   * load parameter where the totla strength is specified not the 
   * energy dependence. THis is useful for OM fits to single energy data
   \param strength is the strength of the surface imaginary potential
   \param R0 is the radius of the potential in fm
   \param a0 is the diffuseness of the potential in fm
  */
void surfaceTF::load(double strength, double R0, double a0)
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

void surfaceTF::SetEnergy(double Ecm0)
{
  Ecm = Ecm0;
  double V = TF.funct(Ecm);


  Imaginary.init(V,R,a);
  V = TF.deltaV(Ecm);
  Dispersive.init(V,R,a);

  V = TF.derDeltaV();
  double VV = TF.derDeltaVParticleHole(); 
  ParticleHole.init(VV,R,a);

  /*
  if (abs(Ecm-Efermi) - 10 < 0.)
    {
      double V1 = TF.deltaV(Efermi+10.);
      double V2 = TF.deltaV(Efermi-10.);
      V = (V1-V2)/20.;
    }
  */
  DerDispersive.init(V,R,a);



  

}
//***************************************************************
/**
 * returns the surface imaginary potential at a given radius.
 * the fubction setEnergy(Ecm) must be run before using this function
\param r is radial distance in fm
 */
double surfaceTF::ImaginaryPot(double r)
{


  return Imaginary.DerWoodSaxon(r);
}
//***************************************************************
  /**
   *returns the real dispersive correction to the surface imaginary
   *potential. The function setEnergy(Ecm) must be run beforehand
   \param r is the radial distance in fm
  */
double surfaceTF::DispersiveCorrection(double r)
{
  return Dispersive.DerWoodSaxon(r);
}
//********************************************************
  /**
   *returns the derivative with respect to energy of the dispersive correction
\param r is the radial distance in fm
   */
double surfaceTF::DerivativeDispersive(double r)
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
double surfaceTF::DerivativeParticleHole(double r)
{
  return ParticleHole.DerWoodSaxon(r);
}
//********************************************
  /**
   * returns the magnitude of the imaginary potential
     \param Ecm in the center-of-mass energy in MeV
   */
double surfaceTF::CentralImaginaryPotential(double Ecm)
{
  return TF.funct(Ecm);
}
//*******************************************
  /**
   * returns the magnitude of the dispersive correction
     \param Ecm in the center-of-mass energy in MeV
   */
double surfaceTF::CentralDeltaRealPotential(double Ecm)
{
  return TF.deltaV(Ecm);
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

double surfaceTF::CentralDerDeltaRealPotential(double Ecm)
{
  double VV = TF.derDeltaV(Ecm);
  /*
  if (abs(Ecm-Efermi)-10. < 0)
    {
      double V1 = TF.deltaV(Efermi+10.);
      double V2 = TF.deltaV(Efermi-10.);
      VV = (V1-V2)/20.; 
    }
  */
  return VV;
}
//****************************************
double surfaceTF::CentralParticleHole(double Ecm)
{

  TF.deltaV(Ecm);
  return TF.derDeltaVParticleHole(); 
}

//***************************************
  /**
   *returns the maximum of the imagimary potential as a function of Ecm
   * units are MeV
  */
double surfaceTF::getMaxW()
{
  double Wmax = 0.;
  for (int i=5;i<50;i++)
    {
      double Ecm = (double)i;
      double W = TF.funct(Ecm);
      if (W > Wmax) Wmax = W;
    }
  return Wmax;
}
