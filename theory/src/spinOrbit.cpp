#include "../include/spinOrbit.h"

/**
 * loads the DOM spin orbit parameters
 */
void spinOrbit::load(double R0, double a0, 
     double Vzero0, double e0, double Efermi0,double AW0, double BW0)
{
  R = R0;
  a = a0;
  Vzero = Vzero0;
  e = e0;
  Efermi = Efermi0;
  AW = AW0;
  BW = BW0;
  Form.init(AW,BW,0.,0.,Efermi,4,0,0.,0.);
  
}
//**********************************************************
/**
 * loads the strength of the real and imaginary spinorbit
 * this version of load is udeful for OM fits to single-energy data
 \param V0 is real strength of spin orbit potential
 \param W0 is imaginary stregth of spin orbit potential
 \param R0 is radius of spin-orbit potential in fm
 \param a0 is the diffuseness of spin-orbit potential in fm
*/
void spinOrbit::load(double V0, double W0, double R0, double a0)
{
  R = R0;
  a = a0;
  V = V0;
  W = W0;
  Real.init(V,R,a);
  Imag.init(W,R,a);
  DerDispersive.init(0.,R,a);
  DerParticleHole.init(0.,R,a);
}

//*******************************************************************
void spinOrbit::SetEnergy(double Ecm0)
{
  Ecm = Ecm0;
  V  = Vzero-e*Ecm+Form.DeltaRealPotential(Ecm);

  Real.init(V,R,a);
  W = Form.ImaginaryPotential(Ecm);
  Imag.init(W,R,a);
  DerV = Form.DerDeltaHole(Ecm) + Form.DerDeltaParticle(Ecm);
  DerDispersive.init(DerV,R,a);


  double DerivParticleHole;
  if (Ecm > Efermi) DerivParticleHole = Form.DerDeltaParticle(Ecm);
  else DerivParticleHole = Form.DerDeltaHole(Ecm);
  DerParticleHole.init(DerivParticleHole,R,a);

}
//***************************************************************
double spinOrbit::RealPotential(double r, double LdotSigma)
{
  return Real.DerWoodSaxon(r)/r/2.*LdotSigma/a;
}

//***************************************************************
double spinOrbit::ImaginaryPotential(double r, double LdotSigma)
{
  return Imag.DerWoodSaxon(r)/r/2.*LdotSigma/a;
}
//***************************************************************
double spinOrbit::DerivativeDispersive(double r, double LdotSigma)
{
  return DerDispersive.DerWoodSaxon(r)/r/2.*LdotSigma/a;
}
//***************************************************************
double spinOrbit::DerivativeParticleHole(double r, double LdotSigma)
{
  return DerParticleHole.DerWoodSaxon(r)/r/2.*LdotSigma/a;
}
