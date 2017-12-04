#include "../include/potPara.h"
#include <iostream>

using namespace std;


//*********************************************************************
// potPara constructor
potPara::potPara( double V0, double R0, double a0)
  :V(V0),R(R0),a(a0)
{}
//********************************************************************
// initialization
void potPara::init(double V0, double R0, double a0)
{
  V = V0;
  R = R0;
  a = a0;
}
//*****************************************************************
// Wood Saxon form factor
double potPara::WoodSaxon(double r)
{
  return -V/(1.+exp((r-R)/a));
}
//******************************************************************
//derivative of Wood Saxon *-4a
double potPara:: DerWoodSaxon(double r)
{
  float fact = exp((r-R)/a);
  return -4.*fact*V/pow(1+fact,2);
}
//******************************************************************
//alternative surface potential - Gaussian
double potPara::Gauss(double r)
{
  return -V*exp(-pow((r-R)/a,2));
}
//*****************************************************************
potPara potPara::copy()
{
  potPara out(V,R,a);
  return out;
}
//******************************************************************
void potPara::Print()
{
  cout << "V= " << V << " R= " << R << " a= " << a << endl;
}

