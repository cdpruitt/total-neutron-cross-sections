#ifndef _whit
#define _whit

#include <iostream>
#include <cmath>
#include <new>
using namespace std;


/**
 *\brief Whittaker functions
 *
 *calculates Whittaker functions - asymptotic for bound-state with 
 *Coulomb interaction
 */

class whit
{
 public:
  whit(int);
  ~whit();

  static double const pi;
  static double const EulerGamma;

  double getWaveFunction(double nu, int l, double rho);
  double AsymptoticExpansion(double,int,double);
  double derivative;  //!< derivative of function
  double ZeroEnergy(double,int);
  double AscendingSeries(double,int,double);
  double whittakerW(double nu,int l,double x);
  double hypergeometricU(double a, double b, double x);
  double gamma2(double x);

  int method; //!< method used to calculate hypergeometric function



 protected:
  void paddie();
  int n;  //!< maximum order
  int m;
  double ** pd;
  double *array;
  double sign(double);
  double oF1(double,double);
  double Kummer (double, double, double, double &);
  double KummerN (int,double,double,double,double &);
  double Pochammer (double, int);
  double psi (double);
  double  gamma (double);

  double chgus(double a, double b, double x);
  double chgul(double a, double b, double x);
  double chgubi(double a, double b, double x);
  double chguit(double a, double b, double x);
  double cchgubi(double a, double b, double x);
  int id; //!< estimated number of significant digits

};

#endif
