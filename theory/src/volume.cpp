#include "../include/volume.h"

/**
 *loads all the parameters needed to calculate potentials in the DOM
 *formalism
 */
void volume::load(double Rzero0, double deltaR0, double expR0, 
                 double a0, double A0,
		 double B0, double Ep0, double Efermi0, 
                 int m0, int Asy0, double alpha0, double Ea0)
{
  Rzero = Rzero0;
  deltaR = deltaR0;
  expR = expR0;
  a = a0;
  A = A0;
  B = B0;
  Ep = Ep0;
  Efermi = Efermi0;
  m = m0;
  Asy = Asy0;
  alpha = alpha0;
  Ea = Ea0;

  volume.init(A,B,0.,Ep,Efermi,m,Asy,alpha,Ea);
  //surface correction to account for change in radius with energy
//if (deltaR > 0.)surface.init(A*deltaR/4./a,B,expR,Ep,Efermi,m,Asy,alpha,Ea,Ep);


  if (deltaR > 0.)
    surface.init(A*deltaR/4./a,B,expR,(double)m,Efermi,Ea,alpha,Ep);
  else surface.init(0.,10.,10.,4.,Efermi,Ea,alpha,Ep);
} 
//**********************************************************
  /**
   * loads simple imaginary potential of const strength, for use in fitting
   * data from one energy in standard OM
   \param strength is the strength of the volume imaginary potential in MeV
   \param R0 is the radius in fm
   \param a0 is the diffuseness in fm
  */
void volume::load(double strength, double R0, double a0)
{
  deltaR = 0.;
  Rzero = R0;
  a = a0;
  ImaginaryVolume.init(strength,Rzero,a0);
  DispersiveVolume.init(0.,Rzero,a);
  DerDispersiveVolume.init(0.,Rzero,a);
  ParticleHoleVolume.init(0.,Rzero,a);
}
//***********************************************************
  /**
   * in the DOM formalism for each given energy this needs to be run
   * before calling the 
   * subsequent functions with give radial dependies
   \param Ecm0 is the center -f mass energy in MeV
   */
void volume::SetEnergy(double Ecm0)
{
  Ecm = Ecm0;


  ImaginaryVolume.init(volume.ImaginaryPotential(Ecm),Rzero,a);

  DispersiveVolume.init(volume.DeltaRealPotential(Ecm),Rzero,a);

  DerDispersiveVolume.init(volume.DerDeltaHole(Ecm)+
                     volume.DerDeltaParticle(Ecm),Rzero,a);


  double R = Rzero;// + deltaR/2.*exp(-pow(B/expR,2)/2.);
  if (deltaR > 0.)
    {
     ImaginarySurface.init(surface.funct(Ecm),R,a);

     DispersiveSurface.init(surface.deltaV(Ecm),R,a);

     DerDispersiveSurface.init(surface.derDeltaV(Ecm),R,a);


    }


  double DerParticleHoleVolume;
  double DerParticleHoleSurface= 0.;
  if (Ecm > Efermi)
    {
     DerParticleHoleVolume = volume.DerDeltaParticle(Ecm);
     if (deltaR > 0.)DerParticleHoleSurface = surface.derDeltaVParticle(Ecm);
    }
  else 
    {
    DerParticleHoleVolume = volume.DerDeltaHole(Ecm);
    if (deltaR > 0.) DerParticleHoleSurface = surface.derDeltaVHole(Ecm);
    }
   ParticleHoleVolume.init(DerParticleHoleVolume,Rzero,a);
  if (deltaR> 0.) ParticleHoleSurface.init(DerParticleHoleSurface,R,a);
}
//***************************************************************
// returns the volume imaginary potential at a given radius
double volume::ImaginaryPot(double r)
{
  double out = ImaginaryVolume.WoodSaxon(r);
  if (deltaR > 0.) out += ImaginarySurface.DerWoodSaxon(r);
  return out;
}
//***************************************************************
//returns the real dispersive correction to the volume imaginary
//potential
double volume::DispersiveCorrection(double r)
{
  double out = DispersiveVolume.WoodSaxon(r);
  if (deltaR > 0.) out  +=  DispersiveSurface.DerWoodSaxon(r);  
  return out;
}
//********************************************************
//returns the derivative with respect to energy of the dispersive correction
double volume::DerivativeDispersive(double r)
{
  double out = DerDispersiveVolume.WoodSaxon(r);
  if (deltaR > 0.) out +=  DerDispersiveSurface.DerWoodSaxon(r);  
  return out;
         
}
//**********************************************************
//returns the derivative of the hole or particle contribution to 
//the dispersive corrections
//for E > Efermi - particle contribution
//for E < Efermi - hole contribution
//used to calculate occupation probabilities
double volume::DerivativeParticleHole(double r)
{
  double out =  ParticleHoleVolume.WoodSaxon(r);
  if (deltaR > 0.) out += ParticleHoleSurface.DerWoodSaxon(r);
   return out;
}
