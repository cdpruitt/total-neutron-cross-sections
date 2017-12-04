#include "disperse.h"


/**
 *!\brief parameterization of the surface potential with gap and 2 fermi functions
 *
 * A parametrization of the energy dependence of the surface imaginary 
 * potential and also gives the dispersive corrections, etc , associated
 * with this
 * if \f$ W(E) = 0 \f$ if \f$|X|-W_{start} < 0 /f$ otherwise
 *\f$ W(E) = A \frac{1}{1+\exp\left( \frac{|X|-C}{D}\right)} \frac{\exp\left(|X|-W_{start}\right)-1}{\exp\left(|X|-W_{start}\right)+1} \f$ ,
 * where \f$ X = E-E_{Fermi} \f$ 
 */
 
class twoFermi : public disperse
{
 protected:
  double A;
  double B;
  double C;
  double D;
  double Ef; //!< Fermi energy in MeV
  double Wstart;

 public:
  twoFermi(){};
  twoFermi(double A, double B, double C, double D, 
            double Wstart,double Ef); 
  void init(double A, double B, double C, double D, 
            double Wstart, double Ef); 

  double funct(double E);
  double derFunct(double E);

  double functX(double x);
  double derFunctX(double x);

  double deltaV(double E);
  double derDeltaV(double E);

  double derDeltaVHole(double E);
  double derDeltaVParticle(double E);

  double derDeltaVParticleHole();
  double derDeltaV();

  double deltaVAboveBelow(double E);

};
