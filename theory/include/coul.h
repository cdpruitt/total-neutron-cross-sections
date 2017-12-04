#ifndef coul_
#define coul_
#include <iostream>
#include <cmath>
#include <complex>
using namespace std;


/**
 *\brief Coulomb Wave functions
 *
 *
 *Coulomb wave functions using continued-fraction evalution of Coulomb
 *Functions and their derivatives
 *Journal of computational Physics 46 (1982) 171
 *The values of F, G and derivative dF, and dG are correct except
 *their sign may be wrong. Either all right or all the wrong sign
 */


class coul
{
 private:
  double BB;
 public:
  double LogDerF(int ,double, double );
  complex<double> LogDerH(int, double, double);
  double abSum (complex<double>);
  int init(int, double, double);
  double F;
  double dF;
  double G;
  double dG;
};
#endif
