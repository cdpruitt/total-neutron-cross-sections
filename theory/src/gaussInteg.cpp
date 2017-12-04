#include "../include/gaussInteg.h"

//**********************************************************************
 gaussInteg::gaussInteg( const gaussInteg &A)
{
  nPoints = A.nPoints;
  pos = A.pos;
  weight = A.weight;
}
//***********************************************************************
 gaussInteg::gaussInteg(int i0, const double *x0,  const double *w0)
{
  nPoints = i0;
  pos = x0;
  weight = w0;
}
//***********************************************************************
//***********************************************************************
//returns the position
double gaussInteg::x(int i) 
{
  return *(pos+i);
}
//***********************************************************************
//returns the weight
double gaussInteg::w(int i) 
{
  return *(weight+i);
}
