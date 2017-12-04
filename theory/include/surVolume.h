#include "disperse.h"

/**
 *\brief parameterization of surface with two Gaussions
 *
 * A parametrization of the energy dependence of the surface imaginary 
 * potential and also gives the dispersive corrections, etc , associated
 * with this
 *\f$ W(E) = \frac{A}{2\sqrt{2\pi}\sigma} \frac{\|X\|^{m}}{\|X\|^{m}+B^m} \left( \exp\left(-\frac{(X-E_0)^2}{2 \sigma^2}\right) + \exp\left(-\frac{(X+E_0)^2}{2 \sigma^2} \right)\right) \f$
 * where \f$ X = E-E_{Fermi} \f$ 
 */

class surVolume : public disperse
{
 protected:
  double A;
  double B;
  double sigma;
  double Ef; //!< Fermi energy in MeV
  double m;
  double Ea;
  double masy;
  double alpha;
  double Ep;
 public:
  surVolume(){};
  surVolume(double A, double B, double sigma, double m, double Ef, 
	    double Ea, double alpha, double Ep); 
  void init(double A, double B, double sigma, double m, double Ef, 
	    double Ea, double alpha, double Ep); 

  double funct(double E);
  double derFunct(double E);

  double functX(double x);
  double derFunctX(double x);

  double deltaV(double E);
  double derDeltaV(double E);

  double derDeltaVHole(double E);
  double derDeltaVParticle(double E);

};
