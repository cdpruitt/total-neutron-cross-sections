#include "../include/asy.h"

double const  asy::pi = acos(-1.);

asy::asy()
{
  Gauss16 = new gauss16();
}
//**************************************************************
asy::~asy()
{
  delete Gauss16;
}
//***********************************************************************
void asy::init(double alpha0, double C0, double Ea0, double Efermi0)
{
  alpha = alpha0;
   C = C0;
   Ea = Ea0;
   Efermi = Efermi0;
   El = Ea + Efermi;
   El32 = pow(El,3./2.);
   El12 = sqrt(El);
}
//**********************************************************
//imaginary potential
//Ex = E - Efermi
double asy::W(double Ex)
{
  double E = Ex + Efermi;
  return alpha*(sqrt(E) + El32/2./E - 3./2.*El12)*exp(-C*fabs(Ex));
}
//***************************************************************
//first derivative of W
double asy::dW(double Ex)
{
  double E = Ex + Efermi;
  double fact= alpha*(1./2./sqrt(E) - El32/2./pow(E,2))*exp(-C*fabs(Ex));
  double plusMinus = 1;
  if (Ex < 0.)  plusMinus = -1.;
  fact += -C*plusMinus*W(Ex);
  return fact;
}
//****************************************************************
//second derivative of W
double asy::d2W(double Ex)
{
  double E = Ex + Efermi;
  double fact = alpha*(-1./4/pow(E,3./2.) + El32/pow(E,3))*exp(-C*fabs(Ex));
  double plusMinus = 1;
  if (Ex < 0.) plusMinus = -1.;
  fact += -2.*C*plusMinus*alpha*(1./2./sqrt(E) - El32/2./pow(E,2))
         *exp(-C*fabs(Ex));
  fact += C*C*W(Ex);
  return fact;
}
//********************************************************************
//dispersive correction
double asy::dispersive(double E)
{
  double Ex = E - Efermi;
  double Ead = fabs(Ex-Ea);

  //limits for integrations 
  double Elarge = Ea;
  double Esmall = Ead;

  if (Ead > Elarge) 
    {
     Elarge = Ead;
     Esmall = Ea;
    }
  double Ehigh;
  double Elow;

  double derivativeOld=0.;
  int Nmax = 0;
  double Emax[3];
  Emax[0] = Elarge;
  for (int i=0;i<15;i++)
    {
      double x = Elarge + (double)i*10.;
      double derivative = (dW(Ex+x) - dW(x))/x  - (W(Ex+x)-W(x))/x/x;
      //cout << x << " " << (W(Ex+x) - W(x))/x << " " << derivative << endl;
      if (i > 0 && derivative*derivativeOld < 0.)
	{

          double xx = x - 5.;
          for (;;)
	    {
             double d =  (dW(Ex+xx) - dW(xx))/xx  - (W(Ex+xx)-W(xx))/xx/xx;
             double d2 = (d2W(Ex+xx)- d2W(xx))/xx - 2.*(dW(Ex+xx)-dW(xx))/xx/xx
	       + 2.*(W(Ex+xx)-W(xx))/pow(xx,3);
             //cout << xx << " " << d << " " << d2 << endl;
             double dx = -d/d2;
             //cout << "dx = " << dx <<endl;
             if (fabs(dx) < .1) break;
	     xx += dx;
	    }
          Nmax++;
	  Emax[Nmax] = xx;
          //cout << Nmax << " " << Emax[Nmax] << endl;
	}
      derivativeOld = derivative; 
    }

  double Integral = 0.;
  double dIntegral = 0.;
  if (Nmax > 0)
    {
      for (int i=0;i<Nmax;i++)
	{
	  Elow = Emax[i];
          Ehigh = Emax[i+1];
          //cout << "Elow= " << Elow << " Ehigh= "<< Ehigh << endl;
          for (int j=0;j<Gauss16->nPoints;j++)
	    {
	      double ee = (Ehigh-Elow)/2.*Gauss16->x(j) + (Ehigh+Elow)/2.;
	      double weight = Gauss16->w(j)*(Ehigh-Elow)/2.;
              double fact = (W(Ex+ee)-W(ee))/ee;
              double dfact = dW(Ex+ee)/ee;
              //cout << ee << " " << fact << " "  << weight << endl;
              Integral += fact*weight;
              dIntegral += dfact*weight;
	    } 
	}
    }

  //now intregrate to infinity
  Ehigh = 1./Emax[Nmax];
  Elow = 0.;
  for (int j=0;j<Gauss16->nPoints;j++)
    {
      double oneOveree = (Ehigh-Elow)/2.*Gauss16->x(j) + (Ehigh+Elow)/2.;
      double weight = Gauss16->w(j)*(Ehigh-Elow)/2.;
      double ee = 1./oneOveree;
      double fact = (W(Ex+ee)-W(ee))/ee;
      double dfact = dW(Ex+ee)/ee;
      //cout << ee << " " << fact << " "  << weight << endl;
      Integral += fact*weight*pow(ee,2);
      dIntegral += dfact*weight*pow(ee,2);
     } 


  //integral from 0 to Esmall
  if (Ex > Ea) 
    {
      Ehigh = Esmall;
      Elow = 0.;
      for (int j=0;j<Gauss16->nPoints;j++)
	{
	  double ee = (Ehigh-Elow)/2.*Gauss16->x(j) + (Ehigh+Elow)/2.;
	  double weight = Gauss16->w(j)*(Ehigh-Elow)/2.;
          double fact = (W(Ex+ee)-W(Ex-ee))/ee;
          double dfact = (dW(Ex+ee)-dW(Ex-ee))/ee;
          //cout << ee << " " << fact << " "  << weight << endl;
          Integral += fact*weight;
          dIntegral += dfact*weight;
        } 
    }


  //integral from Esmall to Elarge
  Ehigh = Elarge;
  Elow = Esmall;
  for (int j=0;j<Gauss16->nPoints;j++)
    {
      double ee = (Ehigh-Elow)/2.*Gauss16->x(j) + (Ehigh+Elow)/2.;
      double weight = Gauss16->w(j)*(Ehigh-Elow)/2.;
      double fact = 0.;
      double dfact = 0.;
      if (ee > Ea - Ex)
	{
          fact +=W(Ex+ee)/ee;
          dfact += dW(Ex+ee)/ee;
	}
      if (ee > Ea) 
	{
         fact -= W(ee)/ee;
         //cout << "yes1 " << -W(ee)/ee << endl; 
	}
      if (ee < Ex - Ea)
	{
        fact -= W(Ex-ee)/ee;
        dfact -= dW(Ex-ee)/ee;
        //cout << "yes2 " << -W(Ex-ee)/ee << endl;
	}
      //cout << ee << " " << fact << " "  << weight << endl;
      Integral += fact*weight;
      dIntegral += dfact*weight;
     } 
  DerDelta = dIntegral/pi;
  return Integral/pi;
}
//***************************************************************
double asy::DerDispersive(double energy)
{
  dispersive(energy);
  return DerDelta;
}
