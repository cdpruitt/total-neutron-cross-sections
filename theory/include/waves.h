#include <iostream>

using namespace std;



/**
 *\brief Coulomb wave functions or Spherical Bessel functions
 *
 *Calculates the regular and irregular solution to the 
 *scattering problem for zero potential. Also give the 
 * derivatives. With the Coulomb force these are the Coulomb
 *wave functions F and G. For neutrons these are Spherical
 * Bessel functions.
 */

class waves
{
 public:
  waves(double,double,int);
  ~waves();
  double  *F; // regular wavefunction
  double *dF; // derivative of regular
  double *G; // irregular wavefunction
  double *dG; // derivative of irregular
  double *Sigma; // Coulomb phase shifts
 private:
  void CoulombWave();
  void SphericalBessel();
  double rho;
  double gamma;
  int lMaximum;
};
//*******************************************************
