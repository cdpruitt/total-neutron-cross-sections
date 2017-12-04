//give the following expoential integrals
//Ei , E1, Si, Ci

#include "../include/expInt.h"

using namespace std;

double const expInt::EulerGamma=0.57721566;

expInt::expInt()
{
}
//***********************************************************************
double expInt::E1(double x)
{
  if (x == 0.) return 1.e300;
  else if (x <= 1.)
    {
      double sum = 1.;
      double r = 1.;
      for (int i=1;i<=25;i++)
	{
	  r *= -(double)i*x/pow((double)i+1.,2);
	  sum += r;
	  if (abs(r) <= abs(sum)*1.e-15) break;
	}
      return -EulerGamma-log(x) + x*sum;
    }
  else 
    {
      int m= 20+(int)(80./x);
      double TO = 0.;
      for (int i=m;i>=1;i--) 
	{
	  double ii = (double)i;
	  TO = ii/(1.+ii/(x+TO));
	}
      TO = 1./(x+TO);
      return exp(-x)*TO;
    }
}
//*********************************************************************
double expInt::Ei(double x)
{
  if (x == 0.) return -1.e300;
  else if (x <= 40.)
    {
      double sum = 1.;
      double r = 1.;
      for (int i=1;i<=100;i++)
	{
	  r *=(double)i*x/pow((double)i+1.,2);
	  sum += r;
	  if (abs(r/sum) <= 1.e-15) break;
	}
      return EulerGamma+log(x)+x*sum;
    }
  else 
    {
      double sum = 1.;
      double r = 1.;
      for (int i=1;i<=20;i++)
	{
	  r *= (double)i/x;
	  sum += r;
	}
      return exp(x)/x*sum;
    }
}
//**************************************************************************
//sine and cosine integrals
void expInt::SiCi(double x)
{
  double const EPS = 1.e-15;
  double x2 = pow(x,2);
  if (x == 0.) 
    {
      Ci = -1e300;
      Si = 0.;
      return;
    }
  else if (x <= 16.)
    {
      double xr = -0.25*x2;
      Ci = EulerGamma + log(x) + xr;
      for (int i=2;i<=40;i++)
	{
	  xr *= -0.5*(double)(i-1)/(double)(i*i*(2*i-1))*x2;
	  Ci += xr;
	  if (abs(xr) < abs(Ci)*EPS) break;
	}
      xr = x;
      Si = x;
      for (int i=1;i<=40;i++)
	{
	  xr *= -0.5*(double)(2*i-1)/(double)(i*(4*i*i+4*i+1))*x2;
	  Si += xr;
	  if (abs(xr) < abs(Si)*EPS) break;
	}
      return;
    }
  else if (x <= 32.)
    {
      int m = (int)(47.2+.82*x);
      double BJ[102];
      double xa1 = 0.;
      double xa0 = 1.e-100;
	for (int i=m;i>=1;i--)
	  {
	    double xa = 4.*(double)i*xa0/x - xa1;
            BJ[i] = xa;
	    xa1 = xa0;
	    xa0 = xa;
	  }
      double xs = BJ[1];
      for (int i=3;i<=m;i+=2) xs +=2.*BJ[i];
      BJ[1] /= xs;
      for (int i=2;i<=m;i++) BJ[i] /= xs;
      double xr = 1.;
      double xg1 = BJ[1];
      for (int i=2;i<=m;i++)
	{
	  xr *= 0.25*pow((double)(2*i-3),2)/(double)((i-1)*(2*i-1)*(2*i-1))*x;
	  xg1 += xr*BJ[i];
	}
      xr = 1.;
      double xg2 = BJ[1];
      for (int i=2;i<=m;i++)
	{
	  xr *= 0.25*pow((double)(2*i-5),2)/(double)((i-1)*(2*i-3)*(2*i-3))*x;
	  xg2 += xr*BJ[i];
	}
      double xcs = cos(x/2.);
      double xss = sin(x/2.);
      Ci = EulerGamma + log(x) - x*xss*xg1+2.*xcs*xg2 - 2.*pow(xcs,2);
      Si = x*xcs*xg1 + 2.*xss*xg2 - sin(x);
    }
  else
    {
      double xr = 1.;
      double xf = 1.;
      for (int i=1;i<=9;i++)
	{
	  xr *= -2.*(double)(i*(2*i-1))/x2;
	  xf += xr;
	}
      xr = 1./x;
      double xg = xr;
      for (int i=1;i<=8;i++)
	{
	  xr *= -2.*(double)(i*(2*i+1))/x2;
	  xg += xr;
	}
      Ci = xf*sin(x)/x - xg*cos(x)/x;
      Si = acos(0.) - xf*cos(x)/x - xg*sin(x)/x;
    }
}
