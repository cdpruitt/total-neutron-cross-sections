#include "../include/standard.h"




standard::standard(double A0, double B0, double C0, double E00, double m0, 
		   double Ef0) : disperse()
{
  init(A0,B0,C0,m0,E00,Ef0);
}
//*****************************************************************
void standard::init(double A0, double B0, double C0,  double E00, double m0, 
double Ef0)
{
  A = A0;
  B = B0;
  C = C0;
  E0 = E00;
  m = m0;
  Ef = Ef0;

}
//*************************************************************
double standard::funct(double E)
{
  double x = E - Ef;
  return functX(x);
}
//*************************************************************
double standard::functX(double x)
{
  double xx = abs(x) - E0;
  if (xx <= 0.) return 0.;
  else return A*pow(xx,m)/(pow(xx,m)+pow(B,m))*exp(-C*xx);
}
//***************************************************************
double standard::derFunct(double E)
{
  double x = E - Ef;
  return derFunctX(x);
}
//***************************************************************
double standard::derFunctX(double x)
{
  double xx = abs(x) - E0;
  if (xx <= 0.) return 0.;  
  double denom = pow(xx,m)+pow(B,m);
  double one = pow(xx,m)/denom;
  double derOne = m/xx*(pow(xx,m)/denom
		  - pow(xx,2*m)/pow(denom,2));
  double two = exp(-C*xx);
  double derTwo = -C*two;
  double out = A*(derOne*two + one*derTwo);
  if (x < 0.) out *= -1.;
  return out;
}
//****************************************************************
double standard::deltaV(double E)
{
  return deltaVX(E-Ef);
}
//*************************************************************
double standard::derDeltaV(double E)
{
  return derDeltaVX(E-Ef);
}
//*************************************************************
double standard::derDeltaVHole(double E)
{
  return derDeltaVHoleX(E-Ef);
}
//*************************************************************
double standard::derDeltaVParticle(double E)
{
  return derDeltaVParticleX(E-Ef);
}
