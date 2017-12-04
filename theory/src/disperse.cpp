#include "../include/disperse.h"
#include <iostream>
using namespace std;

double const disperse::pi=acos(-1.);


/**
 * returns the real dispersive correction to the OM potential (MeV)
\param x is energy relative to Fermi energy in MeV
 */
double disperse::deltaVX(double x)
{


  double sumDispersive = 0.; // dispersive correction
  sumDerivative = 0.; // energy derivative of dispersive correction
  sumPH = 0.; // particle-hole factor
  double yStart;
  if (x > 100.) yStart = x - 100.;
  else yStart = 0.;
  for (int i=0;i<200;i++)
    {
      double y = yStart + (double)i;
      double mag = 0.;
      //corection
      if (x != 0.)
	{
         if (y == 0.) mag = 2.*derFunctX(x);
         else mag = (functX(abs(x+y))-functX(abs(x-y)))/y;

         if (i == 0) mag /= 2.;
	}
      sumDispersive += mag;

      //derivative of correction
      if (y == 0.) y =.001;
      mag = (derFunctX(x+y)-derFunctX(x-y))/y;
      if (i == 0) mag /= 2.;

      sumDerivative += mag;


      // hole or particle factor
      y = (double)i;
      mag = functX(y)/pow(y+abs(x),2);
      if (i == 0) mag /= 2.;
      if (i == 0 && x == 0.) mag = 0.;
      sumPH += mag;


    }
  sumDispersive /= pi;
  sumDerivative /= pi;
  sumPH /= pi;

  return sumDispersive;
}
//*************************************************************************
  /**
   * returns the energy derivative of the dispersive correction
   \param x is energy relative to Fermi energy in MeV
   */
double disperse::derDeltaVX(double x)
{
 
  double sum = 0.;
  double yStart;
  if (x > 100.) yStart = x - 100.;
  else yStart = 0.;
  for (int i=0;i<200;i++)
    {
      double y = yStart + (double)i;
      if (y == 0.) y =.001;
      double mag = (derFunctX(x+y)-derFunctX(x-y))/y;

      if (i == 0) mag /= 2.;

      sum += mag;
    }
  sum /= pi;

  return sum;
 
  return sumDerivative;
}
//************************************************************************
  /**
   * returns the factor used to calculate the hole occupancy in the DOM
   \param x is energy relative to Fermi energy in MeV
   */
double disperse::derDeltaVHoleX(double x)
{
  if (x > 0) abort();
  double sum = 0.;
  for (int i=0;i<200;i++)
    {
      double y = (double)i;
      double mag = functX(y)/pow(y-x,2);
      if (i == 0) mag /= 2.;
      if (i == 0 && x == 0.) mag = 0.;
      sum += mag;
    }
  sum /= pi;
  return sum;
}
//************************************************************************
  /**
   * returns the factor used to calculate the particle occupancy in the DOM
   \param x is energy relative to Fermi energy in MeV
   */
double disperse::derDeltaVParticleX(double x)
{
  return derDeltaVHoleX(-x);
}
//**************************************************
  /**
   * returns the factor used to calculate the hole (particle) occupancy 
   for energies below (above) the Fermi energy. The deltaVX function must
   be run before using this function
   */
double disperse::derDeltaVParticleHoleX()
{
  return sumPH;
}
//***********************************************
  /**
   * returns the derivative of the dispersive correction. 
   * The deltaVX function must be run before using this function 
   */
double disperse::derDeltaVX()
{
  return sumDerivative;
}
//*************************************
  /**
   * if x is greater than zero, returns zero
   * else returns functX
   \param x is energy relative to Fermi energy in MeV
   */
double disperse::functXzeroAbove(double x)
{
  if (x > 0.) return 0.;
  else return functX(x);
}
//*************************************
  /**
   * if x is less than zero, returns zero
   * else returns functX
   \param x is energy relative to Fermi energy in MeV
   */
double disperse::functXzeroBelow(double x)
{
  if (x < 0.) return 0.;
  else return functX(x);
}
//*************************************
  /**
   * if x is greater than zero, returns zero
   * else returns derFunctX
   \param x is energy relative to Fermi energy in MeV
   */
double disperse::derFunctXzeroAbove(double x)
{
  if (x > 0.) return 0.;
  else return derFunctX(x);
}
//*************************************
  /**
   * if x is less than zero, returns zero
   * else returns derFunctX
   \param x is energy relative to Fermi energy in MeV
   */
double disperse::derFunctXzeroBelow(double x)
{
  if (x < 0.) return 0.;
  else return derFunctX(x);
}
/**
 * returns the real dispersive correction to the OM potential (MeV)
\param x is energy relative to Fermi energy in MeV
 */
double disperse::deltaVXAboveBelow(double x)
{


  sumDispersiveAbove = 0.; // dispersive correction from above Ef 
  sumDispersiveBelow = 0.; // dispersive correction from below Ef
  sumDerivativeAbove = 0.; // energy derivative of dispersive correction
  sumDerivativeBelow = 0.; // energy derivative of dispersive correction

  double yStart;
  if (x > 100.) yStart = x - 100.;
  else yStart = 0.;
  for (int i=0;i<200;i++)
    {
      double y = yStart + (double)i;
      double magAbove = 0.;
      double magBelow = 0.;
      //corection
      if (y == 0.) 
        {
	  if (x!= 0)
	    {
             magAbove = 2.*derFunctXzeroBelow(x);
             magBelow = 2.*derFunctXzeroAbove(x);
	    }
        }
      else 
        {
          magAbove = (functXzeroBelow(x+y)-functXzeroBelow(x-y))/y;
          magBelow = (functXzeroAbove(x+y)-functXzeroAbove(x-y))/y;
        }

      if (i == 0) 
        {
          magAbove /= 2.;
          magBelow /= 2.;
        }

      sumDispersiveAbove += magAbove;
      sumDispersiveBelow += magBelow;

      //derivative of correction
      if (y == 0.) y =.001;
      magAbove = (derFunctXzeroBelow(x+y)-derFunctXzeroBelow(x-y))/y;
      magBelow = (derFunctXzeroAbove(x+y)-derFunctXzeroAbove(x-y))/y;
      if (i == 0) 
	{
        magAbove /= 2.;
	magBelow /= 2.;
	}

      sumDerivativeAbove += magAbove;
      sumDerivativeBelow += magBelow;



    }
  sumDispersiveAbove /= pi;
  sumDispersiveBelow /= pi;
  sumDerivativeAbove /= pi;
  sumDerivativeBelow /= pi;

  return sumDispersiveAbove;
}
//**********************************
double disperse::deltaVAbove()
{
  return sumDispersiveAbove - halfIntegral;
}
//**********************************
double disperse::deltaVBelow()
{
  return sumDispersiveBelow + halfIntegral;
}
//**********************************
double disperse::derDeltaVAbove()
{
  return sumDerivativeAbove;
}
//**********************************
double disperse::derDeltaVBelow()
{
  return sumDerivativeBelow;
}
//**********************************
double disperse::findHalfIntegral()
{

  double sum = 0;
  for (int i=1;i<200;i++)
    {
      double y =  (double)i;
      double mag = functX(y)/y;
      sum += mag;
    }
  sum /= 2.;
  halfIntegral = sum;
  return sum;
}
