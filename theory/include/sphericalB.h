#ifndef _sphericalB
#define _sphericalB

#include <iostream>
#include <cmath>

using namespace std;

/**
 *\brief Spherical Bessel functions
 *
 *Calculates spherical Bessal functions
 */

class sphericalB
{

 protected:
  double k0(double);
  double k1(double);
  double k2(double);
  double k3(double);
  double k4(double);
  double k5(double);
  double k6(double);
  double k7(double);

  double y0(double);
  double y1(double);
  double y2(double);
  double y3(double);
  double y4(double);
  double y5(double);
  double y6(double);
  double y7(double);

 public:
  static double const pi;
  double derivative;  //!<derivative of function


  //functions
  double k(int,double);
  double y(int,double);


  double LogDer_k(int,double);
  double LogDer_y(int,double);
};

#endif
