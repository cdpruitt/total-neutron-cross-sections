#include <cmath>
#include <iostream>
#include "../include/legendre.h"

using namespace std;


//this class gives the Legendre Polynomials P0 and P1
//*************************************************

legendre::legendre(int lMaximum0):lMaximum(lMaximum0+1)
{

  //allocate 2d array
  try{
     coef = new double * [lMaximum];
     for (int i=0;i<lMaximum;i++)
       {
         coef[i] = new double [lMaximum];
       }
     }

  catch (bad_alloc &memoryAllocationException)
    {
      cout << "could not allocate Legendre memory "
	   << memoryAllocationException.what() << endl;
    
    }
   //****** zero array
   for (int n=0;n<lMaximum;n++)
     for (int i=0;i<lMaximum;i++) coef[n][i] = 0.;

   //***** load array

   for (int n=0;n<lMaximum;n++)
     {
       int mmMax = n/2;
       for (int mm=0;mm<=mmMax;mm++)
         {
           double sign;
           if (mm%2 == 0) sign=1.;
           else sign = -1.;
              
           int k = n-2*mm;

           //combitorial nCmm
           double comb1 = 1.;
           for (int i=0;i<n-mm;i++) comb1 *= (double)(n-i)/double(n-mm-i);

           //combitorial (2n-2mm)Cn
           double comb2 = 1.;
           for (int i=0;i<(n-2*mm);i++) 
                   comb2 *= (double)(2*n-2*mm-i)/double(n-2*mm-i);    
     
           coef[n][k] = sign*comb1*comb2/pow(2.,n);
	 }
    }


}
//****************************************************************

legendre::~legendre()
{
  for (int i=0;i<lMaximum;i++)
    {
      delete [] coef[i];
    }
  delete [] coef;
}

//****************************************************************
// Legendre polynomial P^0_l(cos(theta))
double legendre::LegendreP0(int l, double theta)
{


  if (l > lMaximum)
    {
      cout << "LegendreP, l > Lmax" << endl;
      return 0;
    }

  double cosTheta = cos(theta);
  double fact = 0.;
  for (int mm = l;mm>=0;mm-=2) fact += coef[l][mm]*pow(cosTheta,mm);
  return fact;
}


//***********************************************************************
//associated Legendre Polynomials P^1_l(cos(theta))
double legendre::LegendreP1(int n, double theta)
{

  if (n <= 0) 
    {
      cout << "associated Legendre polynomial not possible" << endl;
      return 0.;
    }
  if (n > lMaximum)
    { 
     cout << " n> lMaximum in associated Legendre polynomial" << endl;
     return 0.;
    }
  double cosTheta = cos(theta);
  double fact = 0.;
  for (int mm = n;mm>=1;mm-=2) 
     fact += (float)mm*coef[n][mm]*pow(cosTheta,mm-1);
  fact*=-sin(theta);
  return fact; 
}
//***********************************************************************
