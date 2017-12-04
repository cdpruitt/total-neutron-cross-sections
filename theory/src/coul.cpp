#include "../include/coul.h"
using namespace std;




complex<double> coul::LogDerH(int l, double eta, double x)
{
  double const accur = 1e-50;
  int const MaxIterations=2000;
  double const accuh = sqrt(accur);
  complex<double> one(1.,0.);
  complex<double> two(2.,0.);

  complex<double> out(0.,1. - eta/x);

  complex<double> bb(2.*(x-eta),2.);
  complex<double> aa(-pow(eta,2)- pow((double)l,2) - (double)l,eta);
  complex<double> rk(0.,0.);
  complex<double> wi(0.,2.*eta);
  complex<double> rl(0.,1/x);

  if (abSum(bb) < accuh)
    {
      rl *= aa/(aa + rk + wi);
      out +=  rl*(bb+complex<double>(0.,2.));
      aa += 2.*(rk + wi + one);
      bb += complex<double>(0.,4.);
      rk+= complex<double>(4.,0.);
    }
  complex<double> dd = one/bb;
  complex<double> dl = aa*dd*rl;
  for (;;)
    {
      complex<double>outOld = out;
      out += dl;
      rk+= two;
      aa+= rk+wi;
      bb+=complex<double>(0.,2.);
      dd = one/(aa*dd+bb);
      dl*=bb*dd-one;
      //double error = abSum(dl)/abSum(out);
      //cout << out << endl;
      if (rk.real() > 2*MaxIterations) break;
      if (abs(out-outOld) < accur) break;
    }
  out += dl;
  return out;
}

double coul::abSum(complex<double> x)
{
  return abs(x.real()) + abs(x.imag());
}

double coul::LogDerF(int l, double eta, double x)
{
  double const tiny = 1.e-30;
  double const eps = 1.e-30;
  double l1 = (double)l;
  double l2 = l1+1.;
  double S1 = l1/x + eta/l1;
  double S2 = l2/x + eta/l2;
  double BBold = 1.;
  double BBoldOld = 0.;

  double out2 = S2;
  if (out2 == 0.)out2 = tiny;
  double C = out2;
  double D = 0.;
  for (;;)
  {
   double out1 = out2;
   l++;
   l1 = (double)l;
   l2 = l1+1.;
   S1 = S2;
   S2 = l2/x + eta/l2;
   double B = S1 + S2;
   double A = -(1.+pow(eta/l1,2));
   D = B + A*D;
   if (D == 0.) D = tiny;
   C = B + A/C;
   if (C == 0.) C = tiny;
   D = 1./D;
   double delta = C*D;
   out2 = out1*delta;
   BB = B*BBold + A*BBoldOld;
   BBoldOld = BBold;
   BBold = BB;
   //cout << out2 << endl;
   if (abs(delta-1.) < eps)break;
  }
  return out2;
}

int coul::init(int l, double eta, double x)
{
  //double xx = eta + sqrt(pow(eta,2)+ (double)(l*(l+1)));
  //if ( x < xx) cout << " x << xx= " << xx << endl;
  double f = LogDerF(l,eta,x);
  complex<double> h = LogDerH(l,eta,x);
  double p = h.real();
  double q = h.imag();

  //check for JWKB approximation
  if (q < 1.e-10)
    {
      double mu = sqrt((double)(l*(l+1)) + 6./35.);
      double ghalf = sqrt(2.*eta/x + pow(mu/x,2) - 1.);   //ghalf=sqrt(g)
      double S = mu/x + eta/mu;
      double R = sqrt(1.+pow(eta/mu,2));
      double phi2 = 2.*x*ghalf - 2.*eta*atan2(x*ghalf,(x-eta))-
	mu*log(pow((ghalf+S)/R,2)) + log(ghalf);
      double gamma = 2.*ghalf*exp(-phi2);
      double qq = (f-p)/gamma;
      if (qq < 1.e-15) q = qq;
    }

  F = 1./sqrt(pow(f-p,2)/q + q);
  if (BB < 0.) F *= -1.;
  dF = f*F;
  G = (f-p)*F/q;
  dG = p*G - q*F;
  return 1;
}

