#include <cmath>
#include <algorithm>
#include <iostream>
#include "gauss16.h"

using namespace std;

/**
 *\brief energy-asymmetric part of the volume potential
 *
 *The imaginary volume potential is asymmetric in energy about the 
 *Fermi energy. This class defines this asymmetry and gives the 
 *dispersive correction associated with in.
 */

class asy
{
 public:
  double alpha;
  double C;
  double Ea;
  double Efermi;
  double El;
  double El32;
  double El12;
  double DerDelta;

  static double const pi;
 
  asy();
  ~asy();
  double W(double);
  double dW(double);
  double d2W(double); 
  void init(double,double,double,double);
  double dispersive(double);
  double DerDispersive(double);
  gauss16 *Gauss16;
};
