#ifndef _legendre 
#define _legendre
#include <new>


/**
 *\brief Legendre Polynomials P0 and P1
 *
 * calculates the Legendre Polynomials P0 and P1
 */

class legendre
{
public:
  legendre(int);
  ~legendre();
  double LegendreP0(int, double);
  double LegendreP1(int,double);
private:
  int lMaximum;
  double **coef;
};

#endif 
