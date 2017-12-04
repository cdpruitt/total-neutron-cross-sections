#ifndef potpara_
#define potpara_
#include <cmath>


/**
 *\brief Generic radial form factor of potentials
 *
 *Give radial form factors of various type of potentials
 *Wood_saxon, derivative of Wood_saxon, Gaussian
 *
 */

class potPara
{
  public:
  potPara(){};
  potPara(double,double,double);
  void init(double,double,double);  //initialization
  potPara copy();
  double WoodSaxon (double);
  double DerWoodSaxon (double);
  double Gauss (double); // alternative form for surface potential
  void Print();
  double V;
  double R;
  double a;
};
#endif
