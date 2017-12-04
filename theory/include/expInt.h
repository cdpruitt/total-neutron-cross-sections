#ifndef ExpInt_
#define ExpInt_

#include <cmath>
#include <algorithm>


/**
 *\brief exponential Integrals
 *
 *Gives the following expoential integrals Ei , E1, Si, Ci
 */

class expInt
{
 public:
  expInt();
  static double const EulerGamma;
  double Ei(double);
  double E1(double);
  double Si; //Sine Integral
  double Ci; //Cosine Integral
  void SiCi(double);
};
#endif
