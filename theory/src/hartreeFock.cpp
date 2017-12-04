#include "../include/hartreeFock.h"

double const hartreeFock::dd = 5.;
double const hartreeFock::ss = 2.;


/**
 * loads the Hartree Fock class with the energy dependent parameters used in
 * the DOM fits
 */ 
void hartreeFock::load(double R0, double a0, double VHF0, double alpha0,
                        double beta0, double gamma0, double Efermi0,
		       double alphaS0,double betaS0, double gammaS0)
{
  R = R0;
  a = a0;
  VHF = VHF0;
  alpha = alpha0;
  beta = beta0;
  gamma = gamma0;
  Efermi = Efermi0;
  alphaS= alphaS0;
  betaS = betaS0;
  gammaS = gammaS0;
}
//***************************************************************************
  /**
   * loads the Hartree Fock call with fixed energy independent potential
    \param V0 is depth of Hartree Fock potenetial in MeV
    \param dV0 is derivative of depth with respect to com energy
    \param F0 is magnitude of surface real potential in MeV
    \param R0 is radius in fm
    \param a is diffuseness in fm
  */
void hartreeFock::load(double V0, double dV0, double R0, double a0)
{
  R = R0;
  a = a0;
  V = V0;
  dV = dV0;
  RealVolume.init(V,R,a);
  DerivativeRealVolume.init(dV,R,a);
  surface.init(0.,R,a);
  DerivativeSurface.init(0.,R,a);
}
//***************************************************************************
void hartreeFock::SetEnergy(double Ecm)
{
  //  V = VHF - alpha*Ecm + beta*exp(gamma*Ecm);
  //  dV = - alpha + beta*gamma*exp(gamma*Ecm);



  V = VHF - alpha*(Ecm - Efermi)- beta*pow(Ecm-Efermi,2) 
    - gamma*pow(Ecm-Efermi,3);
  dV = -alpha -2.*beta*(Ecm-Efermi) - 3.*gamma*pow(Ecm-Efermi,2);
  RealVolume.init(V,R,a);
  DerivativeRealVolume.init(dV,R,a);
  
  double Vs,dVs;
  double x = Ecm - Efermi - betaS;
  if (x <= 0.) 
    {
      Vs = 0.;
      dVs = 0.;
    }
  else 
    {

      Vs = alphaS*pow(x,2)/(pow(x,2)+pow(gammaS,2));
      dVs = alphaS*2.*pow(gammaS,2)*x/pow(pow(x,2)+pow(gammaS,2),2);
    }


  /*
  if (Ecm-Efermi < -5000.) 
    {
      Vs = 0.;
      dVs = 0.;
    }
  else 
    {
      double x = Ecm - Efermi- gammaS;
      Vs = alphaS*exp(-pow(x/betaS,2)/2.);
      dVs = -Vs*x/pow(betaS,2)/2.;

      // we do not want any of this surface for negative energies
      // use tanh function to suppress this
      double fact = (1.+tanh(Ecm/dd-ss))/2.;
      double dfact = 0.5/dd/pow(cosh(ss-Ecm/dd),2);

      dVs = dVs*fact + Vs*dfact;
      Vs *= fact;



      //Vs = alphaS*pow(Ecm-Efermi,2) + betaS*pow(Ecm-Efermi,3);
      //dVs = 2.*alphaS*(Ecm-Efermi) + 3.*betaS*pow(Ecm-Efermi,2);;
    }
  */
  surface.init(Vs,R,a);
  DerivativeSurface.init(dVs,R,a);
}
//****************************************************************************
double hartreeFock::RealPotential(double r)
{
  return RealVolume.WoodSaxon(r) + surface.DerWoodSaxon(r);
}
//*******************************************************************
double hartreeFock::DerivativePotential(double r)
{
  return DerivativeRealVolume.WoodSaxon(r)+ DerivativeSurface.DerWoodSaxon(r);
}
//*************************************************************
double hartreeFock::RealSurface(double Ecm)
{
     double x = Ecm - Efermi - betaS;
     if (x <= 0. ) return 0.;
     return alphaS*pow(x,2)/(pow(x,2)+pow(gammaS,2));
}
